"""
This module implements the :py:class:`StatusMessage` class for GATHODE/CATHODE.

Growth/CLS Analysis Tool for High-throughput Optical Density
Experiments (GATHODE/CATHODE) StatusMessage class.
"""

# GATHODE  Growth Analysis Tool
#          for High-throughput Optical Density Experiments
#
# Copyright (C) 2014 Nils Christian
#
# This file is part of GATHODE.
#
# GATHODE is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# GATHODE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with GATHODE.  If not, see <http://www.gnu.org/licenses/>.

###########################################
# FIXME this should be inside StatusMessage
def enum(**enums):
    return type('Enum', (), enums)
Severity=enum(message=1,failed=2,warning=3)
#ParameterType=enum(Double=1,UnsignedInteger=2,Integer=3,Boolean=4)
###########################################

class StatusMessage(object):
    """
    Holds (multiple) status information.

    This class allows to collect status information (e.g. for a
    dataset).

    This class can be either a list of statuses or a single status.

    For a single status, the :py:meth:`message
    <.StatusMessage.message>` method returns a string representation
    of the status.

    For a list of statuses, the :py:meth:`message
    <.StatusMessage.message>` method returns a single string,
    consisting of one string representation for each status-type. If a
    status-type exists multiple times, only the highest priority
    message of this type is shown.

    As an alternative to the :py:meth:`message
    <.StatusMessage.message>` method, you can use
    :py:meth:`type2status <.StatusMessage.type2status>` and assemble
    the required string yourself.
    """

    def __init__(self, key=None, shortmsg=None, longmsg=None, severity=Severity.message, **kwargs):
        """
        Initialise with a keyword and a dict of further parameters, or with a list of (sub)statuses.

        :param key: The key-type this single-status belongs to.
        :type key: string
        :param shortmsg: The short message.
        :type shortmsg: string
        :param longmsg: The long message.
        :type longmsg: string
        :param severity: The severity (importance) of this message.
        :type severity: Severity
        :param kwargs: Additional keywords that can be used in shortmsg and longmsg as
                       keyword arguments for `string format
                       <https://docs.python.org/2/library/string.html#string.Formatter.format>`_.
        :type kwargs: string

        """
        if key is None:
            # create an empty list of statuses
            self.substatuses = []
            self.key = None
            self.kwargs_ = None
        elif type(key) is list:
            # create a list of statuses
            self.substatuses = key
            self.key = None
            self.kwargs_ = None
        else:
            self.substatuses = None
            self.key = key
            self.short_ = shortmsg
            self.long_ = longmsg
            self.severity_ = severity
            self.kwargs_ = kwargs

    def _isNone(self):
        """
        Return True if this StatusMessage is None (invalid).

        :return: bool -- True if invalid.

        For internal use only.
        """

        if type(self) is None:
            return True
        if self.key is None and self.substatuses is None:
            return True
        return False

    def isEmpty(self):
        """
        Return True if this is an empty StatusMessage.

        :return: bool -- True if empty.
        """

        if self.substatuses is not None:
            return len(self.substatuses) == 0
        else:
            return self._isNone()

    def statusesWithKey(self, key):
        """
        Return a list of all substatuses that have the given key.

        :param key: The key-type that substatus should match.
        :type key: string

        :return: list(StatusMessage) -- status messages of type key.
        """
        if self.key is not None and self.key == key:
            return [self]
        elif self.substatuses is None:
            return []
        subs=[]
        for s in self.substatuses:
            subs.extend(s.statusesWithKey(key))
        return subs

    def _showRecursive(self,indent=''):
        #
        if type(self) is not None:
            if self.substatuses is not None:
                if len(self.substatuses) == 0:
                    print(indent+str(type(self))+' empty')
                else:
                    print(indent+str(type(self)))
                for s in self.substatuses:
                    if s is not None:
                        s._showRecursive(indent+'    ')
                    else:
                        print(indent+'     None')
            elif self.key is not None:
                print(indent+self.shortmessage())

    def addStatus(self, status):
        """
        Add a substatus to this status.

        :param status: The substatus that should be added.
        :type status: StatusMessage

        Raises exception if this object was not initialised as list of statuses.

        """
        if self.key is not None:
            raise RuntimeError('cannot add status: this is not a list of statuses')
        if status is None:
            self.substatuses.append(status)
        elif type(status) is StatusMessage:
            self.substatuses.append(status)
        elif type(status) is list:
            for s in status:
                self.addStatus(s)
        else:
            raise RuntimeError('neither StatusMessage nor list of StatusMessages: '+str(type(status)))

    def removeStatusesWithKey(self, key):
        """
        Remove all substatuses that have the given key.

        :param key: The key-type that substatus should match.
        :type key: string

        :return: int -- number of removed status messages.
        """
        if self.key is not None and self.key == key:
            self.key = None
            self.short_ = None
            self.long_ = None
            self.severity_ = None
            self.kwargs_ = None
            return 1
        elif self.substatuses is None:
            return 0

        numremoved=0
        newsubstatuses=[]
        for s in self.substatuses:
            subststr=str(s)
            numremoved+=s.removeStatusesWithKey(key)
            if not s._isNone():
                newsubstatuses.append(s)
        self.substatuses=newsubstatuses

        return numremoved

    def shortmessage(self):
        """
        Short message for this status (only works for primary status, not for list of statuses)

        :return: str -- Short message describing this status.
        """
        if self.key is None:
            raise RuntimeError('cannot show shortmessage: this is a list of statuses')
        if self.short_ is not None:
            return self.short_.format(**self.kwargs_)
        return self.key

    def longmessage(self):
        """
        Long message for this status (only works for primary status, not for list of statuses)

        :return: str -- Long message describing this status.
        """
        if self.key is None:
            raise RuntimeError('cannot show longmessage: this is a list of statuses')
        if self.long_ is not None:
            return self.long_.format(**self.kwargs_)
        return self.key

    def messageType(self):
        """
        Message type for this status (only works for primary status, not for list of statuses)

        :return: str -- String describing the type of this status.
        """
        if self.key is None:
            raise RuntimeError('cannot show message type: this is a list of statuses')
        return self.key

    def type2status(self):
        """
        :return: dict -- Dictionary that for each message type contains the highest priority message.
        """
        type2status={}
        if self.key is not None:
            type2status[self.messageType()]=self
            return type2status

        if self.substatuses is None:
            return type2status

        for sts in self.substatuses:
            subtype2status=sts.type2status()
            for msgtype in subtype2status:
                if ((msgtype in type2status)
                    and ((subtype2status[msgtype].severity() <= type2status[msgtype].severity())
                         or (subtype2status[msgtype].severity() == type2status[msgtype].severity()
                             and subtype2status[msgtype].severity2() <= type2status[msgtype].severity2()))):
                    # the priority of the message that already exists is higher, so skip this
                    continue
                type2status[msgtype]=subtype2status[msgtype]
        return type2status

    def message(self,withMsgType=True,seperator='; '):
        """
        :return: str -- Status message.

        For a list of (sub)statuses, only shows the one with the
        highest priority/severity for each category (message
        type).
        """
        msg=''
        type2status=self.type2status()
        # for each type create a single message, the one with the highest priority/severity
        type2statuskeys=list(type2status.keys())
        type2statuskeys.sort()
        for msgtype in type2statuskeys:
            if msg is not '':
                msg+=seperator
            if withMsgType:
                msg+=msgtype+' '
            if type2status[msgtype].severity() is Severity.warning:
                msg+='WARNING '
            msg+=type2status[msgtype].longmessage()
        return msg

    def __str__(self):
        return self.message()

    def severity(self):
        """
        :return: Severity -- The severity/priority of this status.

        If this is a list of (sub)statuses return the higest.
        """
        if self.key is not None:
            return self.severity_
        else:
            sev=Severity.message
            for sts in self.substatuses:
                thissev=sts.severity()
                if thissev > sev:
                    sev=thissev
            return sev

    def severity2(self):
        """
        :return: Severity -- The second-level severity/priority of this status.

        For this to work the parameter dict has to contain the
        'severity2' entry; if this is not set return 0.

        If this is a list of (sub)statuses return the higest.
        """
        if (self.kwargs_ is None) or ('severity2' not in self.kwargs_) or (self.kwargs_['severity2'] is None):
            return 0
        return self.kwargs_['severity2']
