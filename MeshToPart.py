"""Mesh to Part - Fusion 360 add-in. THIS FILE IS THE SHELL, NOT THE WORK.

Lifecycle, the toolbar button, event wiring, and an optional dev bridge. All the
real code lives in mtp_core.py and is RE-IMPORTED FRESH on every command and every
dev job, so editing that file takes effect immediately - no Stop/Run, no restart.
Only changes to this file need the add-in reloaded, so keep it boring.

DEV BRIDGE (off unless the marker file exists)
    Fusion has no headless mode and no CLI, so nothing outside Fusion can run a
    script. Creating devbridge.on next to this file turns on a folder watcher that
    executes mtp_core.DEV_COMMANDS jobs, which is how this add-in gets tested inside
    a live Fusion instead of by hand.

    THREADING RULE: the API is not thread-safe. The poller thread touches only the
    filesystem and hands work to the main thread through a custom event. Never move
    an adsk.* call into _poll().
"""

import adsk.core
import adsk.fusion

import importlib
import json
import os
import sys
import threading
import time
import traceback

CMD_ID = 'MeshToPart_Convert'
CMD_NAME = 'Mesh to Part'
WRAP_CMD_ID = 'MeshToPart_Wrap'
WRAP_CMD_NAME = 'Wrap (offset)'
WRAP_TOOLTIP = (
    'Grow or shrink a body by a distance, through a signed distance field.\n\n'
    'Offset Face and Shell work face by face and give up when an outward offset\n'
    'self-intersects in a concavity. A level set cannot hit that case: the material\n'
    'merges wherever the concave radius is smaller than the offset. The tool then\n'
    'measures the offset it actually achieved.'
)
CMD_TOOLTIP = (
    'Convert a mesh body (imported STL) into a solid or Form part.\n\n'
    'You give a maximum error budget and a maximum gap to ignore; the tool prepares '
    'the mesh, converts it with whatever methods this licence allows, then measures '
    'the real deviation against the original mesh and reports it.'
)

ADDIN_DIR = os.path.dirname(os.path.realpath(__file__))
RES_DIR = os.path.join(ADDIN_DIR, 'resources', 'MeshToPart')
LOG_PATH = os.path.join(ADDIN_DIR, 'MeshToPart.log')

DEV_FLAG = os.path.join(ADDIN_DIR, 'devbridge.on')
DEV_IN = os.path.join(ADDIN_DIR, 'devjobs', 'inbox')
DEV_OUT = os.path.join(ADDIN_DIR, 'devjobs', 'outbox')
DEV_DONE = os.path.join(ADDIN_DIR, 'devjobs', 'done')
DEV_BEAT = os.path.join(ADDIN_DIR, 'devjobs', 'heartbeat.json')
EVENT_ID = 'MeshToPart_DevJob'
POLL_S = 0.6
JOB_TIMEOUT_S = 900.0

_app = None
_ui = None
_handlers = []
_controls = []
_event = None
_job_handler = None
_thread = None
_stop = threading.Event()
_job_done = threading.Event()


def _log(msg):
    try:
        with open(LOG_PATH, 'a', encoding='utf-8') as fh:
            fh.write('%s  [shell] %s\n' % (time.strftime('%H:%M:%S'), msg))
    except Exception:
        pass


def _core():
    """Fresh import of the worker module. This is the hot reload."""
    if ADDIN_DIR not in sys.path:
        sys.path.insert(0, ADDIN_DIR)
    if 'mtp_core' in sys.modules:
        return importlib.reload(sys.modules['mtp_core'])
    return importlib.import_module('mtp_core')


def _fail(where):
    tb = traceback.format_exc()
    _log('%s\n%s' % (where, tb))
    if _ui:
        try:
            core = sys.modules.get('mtp_core')
            brief = core.short_error(tb) if core else tb.strip().splitlines()[-1]
        except Exception:
            brief = tb.strip().splitlines()[-1]
        _ui.messageBox('%s\n\n%s\n\nFull traceback in\n%s' % (where, brief, LOG_PATH),
                       CMD_NAME)


# --------------------------------------------------------------------------- #
# command handlers - thin, they only delegate
# --------------------------------------------------------------------------- #

class _Created(adsk.core.CommandCreatedEventHandler):
    def notify(self, args):
        try:
            cmd = args.command
            _core().on_created(_app, cmd)
            for h, ev in ((_Execute(), cmd.execute),
                          (_Changed(), cmd.inputChanged),
                          (_Validate(), cmd.validateInputs)):
                ev.add(h)
                _handlers.append(h)
        except Exception:
            _fail('Mesh to Part could not open its dialog.')


class _Execute(adsk.core.CommandEventHandler):
    def notify(self, args):
        try:
            _core().on_execute(_app, _ui, args.command.commandInputs)
        except Exception:
            _fail('Mesh to Part failed.')


class _Changed(adsk.core.InputChangedEventHandler):
    def notify(self, args):
        try:
            # args.inputs is the collection the changed input lives in, which for a
            # grouped input is the group's children - always go via the command.
            _core().on_input_changed(_app, args.input.parentCommand.commandInputs,
                                     args.input.id)
        except Exception:
            _log(traceback.format_exc())


class _Validate(adsk.core.ValidateInputsEventHandler):
    def notify(self, args):
        try:
            args.areInputsValid = _core().on_validate(args.inputs)
        except Exception:
            args.areInputsValid = False



class _WrapCreated(adsk.core.CommandCreatedEventHandler):
    def notify(self, args):
        try:
            cmd = args.command
            _core().on_wrap_created(_app, cmd)
            for h, ev in ((_WrapExecute(), cmd.execute),
                          (_WrapChanged(), cmd.inputChanged),
                          (_WrapValidate(), cmd.validateInputs)):
                ev.add(h)
                _handlers.append(h)
        except Exception:
            _fail('Wrap could not open its dialog.')


class _WrapExecute(adsk.core.CommandEventHandler):
    def notify(self, args):
        try:
            _core().on_wrap_execute(_app, _ui, args.command.commandInputs)
        except Exception:
            _fail('Wrap failed.')


class _WrapChanged(adsk.core.InputChangedEventHandler):
    def notify(self, args):
        try:
            _core().on_wrap_input_changed(
                _app, args.input.parentCommand.commandInputs, args.input.id)
        except Exception:
            _log(traceback.format_exc())


class _WrapValidate(adsk.core.ValidateInputsEventHandler):
    def notify(self, args):
        try:
            args.areInputsValid = _core().on_wrap_validate(args.inputs)
        except Exception:
            args.areInputsValid = False


# --------------------------------------------------------------------------- #
# dev bridge
# --------------------------------------------------------------------------- #

class _JobHandler(adsk.core.CustomEventHandler):
    """Runs on the MAIN thread. Every adsk.* call in a job happens here."""

    def notify(self, args):
        job_id = 'unknown'
        job = None
        try:
            job = json.loads(args.additionalInfo)
            job_id = str(job.get('job_id', 'unknown'))
            core = _core()
            fn = core.DEV_COMMANDS.get(job.get('cmd'))
            if fn is None:
                _write_result(job_id, False, error='unknown cmd %r (have: %s)'
                              % (job.get('cmd'), ', '.join(sorted(core.DEV_COMMANDS))))
            else:
                _write_result(job_id, True, result=fn(_app, job))
        except Exception:
            try:
                _write_result(job_id, False, error=traceback.format_exc())
            except Exception:
                pass
        finally:
            _job_done.set()


def _write_result(job_id, ok, result=None, error=None):
    payload = {'job_id': job_id, 'ok': bool(ok), 'ts': time.time(),
               'result': result if result is not None else {}}
    if error:
        payload['error'] = error
    path = os.path.join(DEV_OUT, '%s.json' % job_id)
    tmp = path + '.part'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, indent=2, default=str)
    os.replace(tmp, path)          # atomic: the reader never sees half a file


def _poll():
    """FILESYSTEM ONLY. No adsk.* calls in this thread, ever."""
    while not _stop.is_set():
        try:
            with open(DEV_BEAT, 'w', encoding='utf-8') as fh:
                json.dump({'ts': time.time(), 'pid': os.getpid()}, fh)
        except Exception:
            pass
        try:
            names = sorted(n for n in os.listdir(DEV_IN) if n.endswith('.json'))
        except Exception:
            names = []
        for n in names:
            if _stop.is_set():
                break
            src = os.path.join(DEV_IN, n)
            try:
                with open(src, encoding='utf-8') as fh:
                    job = json.load(fh)
            except Exception:
                try:
                    os.replace(src, os.path.join(DEV_DONE, n + '.unreadable'))
                except Exception:
                    pass
                continue
            job.setdefault('job_id', os.path.splitext(n)[0])
            try:
                os.replace(src, os.path.join(DEV_DONE, n))    # claim before running
            except Exception:
                pass
            _job_done.clear()
            try:
                _app.fireCustomEvent(EVENT_ID, json.dumps(job))
            except Exception:
                continue
            _job_done.wait(JOB_TIMEOUT_S)                     # one job at a time
        _stop.wait(POLL_S)


def _start_bridge():
    global _event, _job_handler, _thread
    for d in (DEV_IN, DEV_OUT, DEV_DONE):
        os.makedirs(d, exist_ok=True)
    try:
        _app.unregisterCustomEvent(EVENT_ID)      # survive a reload
    except Exception:
        pass
    _event = _app.registerCustomEvent(EVENT_ID)
    _job_handler = _JobHandler()
    _event.add(_job_handler)
    adsk.autoTerminate(False)
    _stop.clear()
    _thread = threading.Thread(target=_poll, daemon=True)
    _thread.start()
    _log('dev bridge watching %s' % DEV_IN)


def _stop_bridge():
    global _event, _job_handler, _thread
    _stop.set()
    _job_done.set()
    if _thread:
        _thread.join(timeout=3.0)
        _thread = None
    try:
        if _event and _job_handler:
            _event.remove(_job_handler)
        _app.unregisterCustomEvent(EVENT_ID)
    except Exception:
        pass
    _event = None
    _job_handler = None


# --------------------------------------------------------------------------- #
# add-in entry points
# --------------------------------------------------------------------------- #

def run(context):
    global _app, _ui
    try:
        _app = adsk.core.Application.get()
        _ui = _app.userInterface

        old = _ui.commandDefinitions.itemById(CMD_ID)
        if old:
            old.deleteMe()
        cmd_def = _ui.commandDefinitions.addButtonDefinition(
            CMD_ID, CMD_NAME, CMD_TOOLTIP, RES_DIR)
        h = _Created()
        cmd_def.commandCreated.add(h)
        _handlers.append(h)

        # the panel list lives in the core so it can be corrected without a reload
        placed = _core().place_button(_ui, cmd_def)
        _log('button placed in %s' % (placed or 'NO PANEL'))

        old_wrap = _ui.commandDefinitions.itemById(WRAP_CMD_ID)
        if old_wrap:
            old_wrap.deleteMe()
        wrap_def = _ui.commandDefinitions.addButtonDefinition(
            WRAP_CMD_ID, WRAP_CMD_NAME, WRAP_TOOLTIP, RES_DIR)
        wh = _WrapCreated()
        wrap_def.commandCreated.add(wh)
        _handlers.append(wh)
        wplaced = _core().place_button(_ui, wrap_def,
                                       _core().WRAP_PANEL_GROUPS)
        _log('wrap button placed in %s' % (wplaced or 'NO PANEL'))

        if os.path.exists(DEV_FLAG):
            _start_bridge()
        _log('started')
    except Exception:
        _fail('Mesh to Part failed to start.')


def stop(context):
    try:
        _stop_bridge()
        for ctrl in _controls:
            try:
                ctrl.deleteMe()
            except Exception:
                pass
        _controls[:] = []
        for cid in (CMD_ID, WRAP_CMD_ID):
            cmd_def = _ui.commandDefinitions.itemById(cid)
            if cmd_def:
                cmd_def.deleteMe()
        _handlers[:] = []
        _log('stopped')
    except Exception:
        _fail('Mesh to Part failed to stop.')
