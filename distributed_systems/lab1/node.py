import sys
import pika
import datetime
import time
from threading import Thread, Lock

def Log(idx, message):
    print("{} [id = {}] {}".format(datetime.datetime.now(), idx, message))

def LogMethod(idx, typename, method):
    Log(idx, "[{}] {}".format(typename, method))

class TMsg:
    def __init__(self, msg_type = "NONE", sender = -1, receiver = -1, numbers = []):
        self.msg_type = msg_type
        self.sender = sender
        self.receiver = receiver
        self.numbers = numbers

    def ToString(self):
        serialized = "{} {} {}".format(self.msg_type, self.sender, self.receiver)
        for x in self.numbers:
            serialized += " {}".format(x)
        return serialized

    def FromString(self, serialized):
        data = serialized.split(' ')
        self.msg_type = data[0]
        self.sender = int(data[1])
        self.receiver = int(data[2])
        if len(data) >= 4:
            self.numbers = [int(x) for x in data[3:]]
        else:
            self.numbers = []


class TRabbitClient:
    def __init__(self, idx):
        self.idx = idx
        self.connection = pika.BlockingConnection(pika.ConnectionParameters('localhost'))
        self.input = self.connection.channel()
        self.output = self.connection.channel()

    def QueueName(self, to_id):
        return "queue_to_{}".format(to_id)

    def SendMessage(self, message, to_id):
        #LogMethod(self.idx, "TRabbitClient", "SendMessage")
        target = self.QueueName(to_id)
        self.output.queue_declare(queue=target)
        self.output.basic_publish(exchange='', routing_key=target, body=message)

    def GetNextMessageForId(self, recipient_id, timeout_seconds = 5):
        #LogMethod(self.idx, "TRabbitClient", "GetNextMessageForId")
        target = self.QueueName(recipient_id)
        self.input.queue_declare(queue=target)
        frame, properties, body = next(self.input.consume(queue=target, inactivity_timeout=timeout_seconds))
        if frame and properties and body:
            self.input.basic_ack(frame.delivery_tag)
            return body.decode()
        LogMethod(self.idx, "ERROR", "LEADER FAILED")
        return None

class TProcessor:
    def __init__(self, ring_size, idx, leader):
        self.mutex = Lock()
        self.ring_size = ring_size
        self.idx = idx
        self.leader = leader
        self.nxt = self.GetNextNode(self.idx)
        self.prev = -1

    def Ping(self):
        LogMethod(self.idx, "TProcessor", "Ping")
        self.mutex.acquire()
        msg = TMsg("PING", self.idx, self.nxt, [])
        self.mutex.release()
        return msg

    def RecvNext(self, serialized):
        LogMethod(self.idx, "TProcessor", "RecvNext Message: {}".format(serialized))
        if serialized:
            return self.ProcessMessage(serialized)
        else:
            if self.prev == self.leader or self.prev == -1:
                return self.InitiateElection()
        return None

    def ProcessMessage(self, serialized):
        LogMethod(self.idx, "TProcessor", "ProcessMessage")
        message = TMsg()
        message.FromString(serialized)
        if message.msg_type == "ELECTION":
            return self.ProcessElection(message)
        if message.msg_type == "LEADER":
            return self.UpdateLeader(message)
        if message.msg_type == "PING":
            return self.ProcessPing(message)
        raise Exception("[{}] Unknown message type! Serialized message: \"{}\"".format(idx, serialized))

    def GetNextNode(self, cur):
        LogMethod(self.idx, "TProcessor", "GetNextNode")
        return (cur + 1) % self.ring_size

    def InitiateElection(self):
        LogMethod(self.idx, "TProcessor", "InitiateElection")
        LogMethod(self.idx, "ELECTION", "START")
        return TMsg("ELECTION", self.idx, self.nxt, [self.idx])

    def ProcessElection(self, message):
        LogMethod(self.idx, "TProcessor", "ProcessElection")
        if self.leader == self.nxt and self.leader != self.idx:
            self.mutex.acquire()
            self.nxt = self.GetNextNode(self.nxt)
            self.mutex.release()
        if self.idx in message.numbers:
            return self.SelectLeader(message)
        else:
            return self.AddSelfAndForward(message)

    def AddSelfAndForward(self, message):
        LogMethod(self.idx, "TProcessor", "AddSeldAndForward")
        message.receiver = self.nxt
        message.numbers.append(self.idx)
        LogMethod(self.idx, "ELECTION", "ADD SELF AND FORWARD MESSAGE")
        return message

    def SelectLeader(self, message):
        LogMethod(self.idx, "TProcessor", "SelectLeader")
        self.leader = max(message.numbers)
        LogMethod(self.idx, "ELECTION", "NEW LEADER {}".format(self.leader))
        return TMsg("LEADER", self.idx, self.nxt, [self.idx, self.leader])

    def UpdateLeader(self, message):
        LogMethod(self.idx, "TProcessor", "UpdateLeader")
        initiator, new_leader = message.numbers
        if new_leader != self.leader:
            self.leader = new_leader
            LogMethod(self.idx, "ELECTION", "NEW LEADER {}".format(self.leader))
        if initiator != self.idx:
            message.receiver = self.nxt
            return message
        LogMethod(self.idx, "ELECTION", "END")
        return None

    def ProcessPing(self, message):
        LogMethod(self.idx, "TProcessor", "ProcessPing")
        self.prev = message.sender
        return None

class TReader(Thread):
    def __init__(self, node):
        Thread.__init__(self)
        self.client = TRabbitClient(node.idx)
        self.node = node

    def run(self):
        while True:
            if self.node.idx != 9:
                serialized = self.client.GetNextMessageForId(self.node.idx)
                response = self.node.RecvNext(serialized)
                if response:
                    self.client.SendMessage(response.ToString(), response.receiver)

class TPinger(Thread):
    def __init__(self, node):
        Thread.__init__(self)
        self.client = TRabbitClient(node.idx)
        self.node = node

    def run(self):
        while True:
            if self.node.idx != 9:
                message = self.node.Ping()
                self.client.SendMessage(message.ToString(), message.receiver)
            time.sleep(1)

class TNode:
    def __init__(self, ring_size, idx, leader):
        self.idx = idx
        self.node = TProcessor(ring_size = ring_size, idx = idx, leader = leader)

        self.pinger = TPinger(node = self.node)
        self.reader = TReader(node = self.node)

        Log(self.idx, "Initialized!")

    def start(self):
        Log(self.idx, "Started!")
        self.pinger.start()
        self.reader.start()

if __name__ == "__main__":
    if (len(sys.argv) < 3):
        print("Use: python3 node.py ring_size leader")
        exit(0)

    ring_size = int(sys.argv[1])
    leader = int(sys.argv[2])

    nodes = [TNode(ring_size = ring_size, idx = i, leader = leader) for i in range(ring_size)]

    for node in nodes:
        node.start()
