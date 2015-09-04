##Brokersettings.
BROKER_URL = 'amqp://myuser:mypassword@localhost:20100/myvhost'
#BROKER_URL = 'amqp://:guest@localhost:20100//'


## Using the database to store task state and results.
CELERY_RESULT_BACKEND = 'amqp://myuser:mypassword@localhost:20110/myvhost'

CELERY_ANNOTATIONS = {'tasks.add': {'rate_limit': '10/s'}}
