# syntax=docker/dockerfile:1
   
FROM node:18-alpine
RUN apk add python3 py3-pip
# RUN pip install pairpro
COPY ./pairpro/ /pairpro/
COPY ./scripts/ /scripts_in_image/
# change work dir heere: WORKDIR /pairpr 
CMD ["python", "/scripts_in_image/train_model.py"]
EXPOSE 3000