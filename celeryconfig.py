##Brokersettings.
BROKER_URL = 'amqp://myuser:mypassword@localhost:20100/myvhost3'#rabbitMQ can handle multiple programs running through it at the same time, 
                                                                #but they all need their own virutal host, if you want to run multiple gas on the same comp these need to be set up
#BROKER_URL = 'amqp://:guest@localhost:20100//'


## Using the database to store task state and results.
CELERY_RESULT_BACKEND = 'amqp://myuser:mypassword@localhost:20110/myvhost3'#again making sure the vhosts are the same is important otherwise results can get mixed up

CELERY_ANNOTATIONS = {'tasks.add': {'rate_limit': '10/s'}}#this controls how fast tasks are moved into the alloted queues
