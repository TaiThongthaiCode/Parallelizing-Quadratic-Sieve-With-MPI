from flask import Flask
from flask_restful import Resource, Api
from math import *

app = Flask(__name__)
api = Api(app)

class FBS(Resource):
    def get(self, n):

        n = int(n)
        n = pow(exp(sqrt(log(n)*log(log(n)))),sqrt(2))
        n = ceil(n)

        return str(n)


api.add_resource(FBS, '/<string:n>')

if __name__ == '__main__':
    app.run()
