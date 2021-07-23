FROM ubuntu:latest

RUN apt-get update && apt-get install make build-essential git python3 -y

RUN useradd -ms /bin/bash bufet2

USER bufet2

WORKDIR /home/bufet2/

COPY . /home/bufet2

# RUN mv /home/bufet2/BUFET2/* /home/bufet2/

# RUN rm BUFET2 -rf

RUN make



CMD /usr/bin/python3 bufet2.py