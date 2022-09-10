import sys
import logging
from os.path import basename

#from gemini_calmgr.utils.dbtools import REQUIRED_TAG_DICT
from recipe_system.cal_service.localmanager import extra_descript, args_for_cals, LocalManager
from recipe_system.cal_service.calrequestlib import get_cal_requests
from gemini_calmgr.cal import get_cal_object

from sqlalchemy.sql.elements import BooleanClauseList, BinaryExpression, Grouping

import astrodata
import ghostdr

from gemini_calmgr.orm.header import Header
from gemini_calmgr.orm.diskfile import DiskFile
from gemini_calmgr.orm.ghost import Ghost

import ghostdr.ghost.primitives_ghost

import re

import sys
import logging
from os.path import basename

from gemini_calmgr.orm import sessionfactory
from recipe_system.cal_service.localmanager import extra_descript, args_for_cals
from gemini_calmgr.cal import get_cal_object


def _check_equals_true(val, calval):
    if calval is True:
        return "pass"
    else:
        return "fail"


def _check_equals_false(val, calval):
    if calval is False:
        return "pass"
    else:
        return "fail"


def _check_equals(val, calval):
    if calval == val:
        return "pass"
    else:
        return "fail"


def _check_not_equals(val, calval):
    if calval != val:
        return "pass"
    else:
        return "fail"


def _check_greater_than(val, calval):
    if calval is not None and val < calval:
        return "pass"
    else:
        return "fail"


def _check_less_than(val, calval):
    if calval is not None and val > calval:
        return "pass"
    else:
        return "fail"


def _check_greater_than_or_equal(val, calval):
    if calval is not None and val <= calval:
        return "pass"
    else:
        return "fail"


def _check_less_than_or_equal(val, calval):
    if calval is not None and val >= calval:
        return "pass"
    else:
        return "fail"


def _check_like(val, calval):
    if len(val) > 2 and '%' in val[1:-1]:
        return "unknown"
    if val is None:
        return "fail"
    if val == '':
        if calval == '':
            return "pass"
        else:
            return "Fail"
    if val == '%' or val == '%%':
        return "pass"
    if val.startswith('%') and val.endswith('%'):
        if calval is not None and val[1:-1] in calval:
            return "pass"
        else:
            return "fail"
    elif val.startswith('%'):
        if calval is not None and calval.endswith(val[1:]):
            return "pass"
        else:
            return "fail"
    elif val.endswith('%'):
        if calval is not None and calval.startswith(val[:-1]):
            return "pass"
        else:
            return "fail"
    if val == calval:
        return "pass"
    else:
        return "fail"


def _check_is_null(val, calval):
    if calval is None:
        return "pass"
    else:
        return "fail"


_re_equals_true = re.compile(r'\w+ = true')
_re_equals_false = re.compile(r'\w+ = false')
_re_equals = re.compile(r'\w+ = :\w+')
_re_not_equals = re.compile(r'\w+ != :\w+')
_re_greater_than = re.compile(r'\w+ > :\w+')
_re_less_than = re.compile(r'\w+ < :\w+')
_re_greater_than_or_equal = re.compile(r'\w+ >= :\w+')
_re_less_than_or_equal = re.compile(r'\w+ <= :\w+')
_re_like = re.compile(r'\w+ LIKE :\w+')
_re_is_null = re.compile(r'\w+ IS NULL')

_checks = [
    (_re_equals_false, _check_equals_false),
    (_re_equals_true, _check_equals_true),
    (_re_equals, _check_equals),
    (_re_not_equals, _check_not_equals),
    (_re_greater_than, _check_greater_than),
    (_re_less_than, _check_less_than),
    (_re_greater_than_or_equal, _check_greater_than_or_equal),
    (_re_less_than_or_equal, _check_less_than_or_equal),
    (_re_like, _check_like),
    (_re_is_null, _check_is_null)
]


def get_status(val, calval, expr):
    for check in _checks:
        if check[0].search(expr):
            return check[1](val, calval)
    return "unknown"


def get_calibration_type(obj):
    if isinstance(obj, Header):
        observation_type = obj.observation_type
        types = obj.types
    else:
        observation_type = obj.observation_type()
        types = obj.tags

    def add_processed(retval, types):
        if 'PROCESSED' in types:
            return True, retval
        return False, retval

    if observation_type == 'FLAT':
        return add_processed('flat', types)
    if observation_type == 'ARC':
        return add_processed('arc', types)
    if observation_type == 'BIAS':
        return add_processed('bias', types)
    if observation_type == 'DARK':
        return add_processed('dark', types)
    if observation_type == 'STANDARD':
        return add_processed('standard', types)
    if 'SLITILLUM' in types:
        return add_processed('slitillum', types)
    return None



def show_line(table_name, key, cal_value, value, expr):
    ascii_codes = {
        "pass": ("\u001b[32m", "\u001b[37m"),
        "fail": ("\u001b[31m", "\u001b[37m"),
        "unknown": ("\u001b[33m", "\u001b[37m")
    }
    status = get_status(value, cal_value, expr)
    start_code, stop_code = ascii_codes[status]
    if (not isinstance(cal_value, str) or len(cal_value) <= 28) \
            and (not isinstance(value, str) or len(value) <= 28):
        print("%s%9s | %18s | %30s | %30s | %s%s" % (start_code, table_name, key, cal_value, value, expr, stop_code))
    else:
        print("%s%9s | %18s | cal: %58s | %s%s" % (start_code, table_name, key, cal_value, expr, stop_code))
        print("%s%9s | %18s | val: %58s | %s%s" % (start_code, '', '', value, '', stop_code))


def debug_binary_expression(clause, cal_obj, header, diskfile, instr):
    if hasattr(clause.left, 'table'):  # isinstance(clause.left, AnnotatedColumn):
        table = clause.left.table
        key = clause.left.key
        val = clause.right.value if hasattr(clause.right, 'value') else None
        if val is None:
            if hasattr(clause.right, 'clauses') and len(clause.right.clauses) > 0:
                vals = []
                for cl in clause.right.clauses:
                    if hasattr(cl, 'value') and cl.value is not None:
                        vals.append("%s" % cl.value)
                val = ', '.join(vals)
            else:
                val = ''
        expr = "%s" % clause
        if table.name == 'header':
            show_line(table.name, key, getattr(header, key), val, expr)
        elif table.name == 'diskfile':
            show_line(table.name, key, getattr(diskfile, key), val, expr)
        else:
            show_line(table.name, key, getattr(instr, key), val, expr)


def debug_boolean_clause_list(clause, cal_obj, header, diskfile, instr, is_or=False):
    for clause in clause.clauses:
        debug_dispatch(clause, cal_obj, header, diskfile, instr)
        # yield x


def debug_dispatch(clause, cal_obj, header, diskfile, instr):
    if isinstance(clause, BooleanClauseList):
        debug_boolean_clause_list(clause, cal_obj, header, diskfile, instr)
    elif isinstance(clause, BinaryExpression):
        debug_binary_expression(clause, cal_obj, header, diskfile, instr)
    elif isinstance(clause, Grouping):
        if 'OR' in str(clause) and isinstance(clause.element, BooleanClauseList):
            # ew, need to debug an OR
            print("\u001b[33mOR Expression:\u001b[37m")
            debug_boolean_clause_list(clause.element, cal_obj, header, diskfile, instr, is_or=True)
            print("\u001b[33mOR Expression Complete\u001b[37m")
        elif 'LIKE' in str(clause):
            debug_binary_expression(clause, cal_obj, header, diskfile, instr)
        else:
            print("Unsupported query element: %s" % str(clause))


def debug_parser(query, cal_obj, header, diskfile, instr):
    for clause in query.query.whereclause.clauses:
        debug_dispatch(clause, cal_obj, header, diskfile, instr)


def build_descripts(rq):
    descripts = rq.descriptors
    for (type_, desc) in list(extra_descript.items()):
        descripts[desc] = type_ in rq.tags
    return descripts


# TODO maybe make the dbtools behavior configurable to allow this without duplication
# This is a hack that convinces the cal code to allow ingests of unprocessed calibrations for the DB
# useful for testing matches
#REQUIRED_TAG_DICT["__dummy__"] = []


def why_not_matching(filename, processed, cal_type, calibration):
    try:
        filead = astrodata.open(filename)
    except Exception as ex:
        logging.error(f"Unable to open {filename} with DRAGONS")
        exit(1)
    try:
        calad = astrodata.open(calibration)
        if cal_type == "auto":
            processed, cal_type = get_calibration_type(calad)
    except:
        logging.error(f"Unable to open {calibration} with DRAGONS")
        exit(2)
    try:
        mgr = LocalManager(":memory:")
        mgr.init_database(wipe=True)
    except:
        logging.error("Unable to setup in-memory calibration manager")
        exit(3)
    try:
        mgr.ingest_file(calibration)
    except Exception as ingestex:
        logging.error("Unable to ingest calibration file")
        raise
        exit(4)

    rqs = get_cal_requests([filead,], cal_type, procmode=None)
    if not rqs:
        logging.error("Unexpected error creating cal requests")
        exit(5)

    reasons = list()
    for idx in range(len(rqs)):
        rq = rqs[idx]
        descripts = build_descripts(rq)
        types = rq.tags
        cal_obj = get_cal_object(mgr.session, filename=None, header=None,
                                 descriptors=descripts, types=types, procmode=rq.procmode)
        method, args = args_for_cals.get(cal_type, (cal_type, {}))

        # Obtain a list of calibrations and check if we matched
        args["return_query"] = True
        if processed:
            args["processed"] = True

        if not hasattr(cal_obj, method):
            print(f"Instrument {calad.instrument()} has no matching rule for {cal_type}")
        else:
            cals, query_result = getattr(cal_obj, method)(**args)

            for cal in cals:
                if cal.diskfile.filename == basename(calibration):
                    logging.info("Calibration matched")
                    exit(0)

            header = mgr.session.query(Header).first()
            diskfile = mgr.session.query(DiskFile).first()

            if calad.instrument().lower().startswith("gmos"):
                instr = mgr.session.query(Gmos).first()
            elif calad.instrument().lower() == "f2":
                instr = mgr.session.query(F2).first()
            elif calad.instrument().lower() == "nifs":
                instr = mgr.session.query(Nifs).first()
            elif calad.instrument().lower() == "niri":
                instr = mgr.session.query(Niri).first()
            elif calad.instrument().lower() == "gnirs":
                instr = mgr.session.query(Gnirs).first()
            elif calad.instrument().lower() == "ghost":
                instr = mgr.session.query(Ghost).first()
            elif calad.instrument().lower() == "nici":
                instr = mgr.session.query(Nici).first()
            elif calad.instrument().lower() == "michelle":
                instr = mgr.session.query(Michelle).first()
            elif calad.instrument().lower() == "gsaoi":
                instr = mgr.session.query(Gsaoi).first()
            elif calad.instrument().lower() == "gpi":
                instr = mgr.session.query(Gpi).first()

            print('Relevant fields from calibration:\n')
            print('Table     | Key                | Cal Value                      '
                  '| Value                          | Expr')
            print('----------+--------------------+--------------------------------'
                  '+--------------------------------+-------------------')
            debug_parser(query_result, cal_obj, header, diskfile, instr)

    if reasons:
        logging.info(reasons)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        logging.error("Useage: why_not_matching <filename> <cal_type> <calibrationfilename>")
    filename = sys.argv[1]
    cal_type = sys.argv[2]
    if cal_type.startswith('processed_'):
        processed = True
        cal_type = cal_type[10:]
    else:
        processed = False
    calibration = sys.argv[3]

    why_not_matching(filename, processed, cal_type, calibration)

