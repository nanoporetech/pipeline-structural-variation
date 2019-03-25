from __future__ import print_function
import json
import logging

class Record(object):

    def __init__(self):
        self.record = dict()

        self.record['retvalue'] = None
        self.record['exit_status'] = None

    def set(self, key, value):
        self.record[key] = value

    def __repr__(self):
        return json.dumps(self.record)


class Telemetry(object):

    EXCEPTION_THROWN = 'Exception thrown'
    WORKFLOW_SUCCESSFUL = 'Workflow successful'

    PASS = "PASS"
    FAIL = "FAIL"

    def __init__(self, filename=None, append=False):
        self.filename = "{}.data.json".format(filename)
        # This will create it if it does not exist already
        if append:
            self.out = open(self.filename, 'a')
        else:
            if not filename:
                self.out = open('/dev/stdout', 'w')
            else:
                self.out = open(self.filename, 'w')

    def __enter__(self):
        self.record = Record()
        return self.record

    def __exit__(self, type, value, traceback):
        if type:
            self.record.set('retvalue', self.EXCEPTION_THROWN)
            self.record.set('exit_status', self.FAIL)

            exception = dict()
            exception['type'] = type.__name__
            exception['message'] = value.message

            self.record.set('exception', exception)

            logging.error(value.message)
        else:
            self.record.set('retvalue', self.WORKFLOW_SUCCESSFUL)
            self.record.set('exit_status', self.PASS)

        print(self.record, file=self.out)

    def __del__(self):
        self.out.close()

