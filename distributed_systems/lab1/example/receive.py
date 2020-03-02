import pika

if __name__ == "__main__":
    connection = pika.BlockingConnection(pika.ConnectionParameters(host='localhost'))
    channel = connection.channel()
    channel.queue_declare(queue='hello')
    print("waiting for messages...")
    def callback(ch, method, properties, body):
        print("received message {}".format(body))
    channel.basic_consume(on_message_callback=callback, queue='hello', auto_ack=True)
    channel.start_consuming()

