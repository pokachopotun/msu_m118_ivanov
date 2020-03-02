import sys
import pika

class TMsg:
    def ToString(self):
        serialized = "{} {}".format(self.msg_type, self.sender)
        for x in self.numbers:
            serialized += str(x)
        return serialized

    def FromString(self, serialized):
        data = serialized.split(' ')
        self.msg_type = data[0]
        self.sender = int(data[1])
        self.numbers = data[2:]

def queue_name(from_id, to_id):
    return "queue_from_{}_to_{}".format(from_id, to_id)

class TNode:
    def __init__(self, ring_size, idx, leader):
        self.ring_size = ring_size
        self.idx = idx
        self.leader = leader
        self.nxt = GetNextNode(self.idx)
        self.last_leader_message = -1

    def SendNext(self, message):
        message.sender = self.idx
        # TODO:
        pass

    def RecvNext(self):
        # TODO:
        pass

    def ProcessMessage(self, serialized):
        msg = TMsg()
        message.FromString(serialized)
        if message.msg_type == "ELECTION":
            ProcessElection(message)
        if message.msg_type == "LEADER":
            ProcessLeader(message)
        if message.msg_type == "DISCOVER":
            ProcessDiscover(message)
        raise Exception("[{}] Unknown message type! Serialized message: \"{}\"".format(idx, serialized))

    def GetNextNode(self, cur):
        return (cur + 1) % self.ring_size

    def ProcessElection(self, message):
        if leader == self.nxt:
            self.nxt = GetNextNode(self.nxt)
        if idx in message.numbers:
            leader = max(message.numbers)
            LeaderSelected(leader)
        else:
            message.numbers.append(self.idx)
            SendNext(message)

    def LeaderSelected(self, leader):
        msg = TMsg()
        msg.msg_type = "LEADER"
        msg.numbers = [self.idx, leader]
        SendNext(msg)

    def ProcessLeader(self, message):
        initiator, self.leader = message.numbers
        if initiator != self.idx:
            SendNext(message)

    def ProcessDiscover(self, message):
        if idx in message.numbers:
            leader = max(message.numbers)
            LeaderSelected(leader)
        else:
            message.numbers.append(self.idx)
            SendNext(message)


if __name__ == "__main__":
    if (len(sys.argv) < 1):
        print("Use python3 node.py idx")
    idx = int(sys.argv[1])
    print("Started node index {}".format(idx))
    print(queue_name(1,2))
#    connection = pika.BlockingConnection(pika.ConnectionParameters('localhost'))
#    channel = connection.channel()
#    channel.queue_declare(queue='hello')
#    channel.basic_publish(exchange='', routing_key='hello', body='hello_world')
#    print("message sent")
#    connection.close()
