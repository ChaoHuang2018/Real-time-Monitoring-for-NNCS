make && \
    ./benchmark1_relu 0.1 50 4 6 1 && \
    ./benchmark1_relu_tanh 0.1 35 4 6 1 && \
    ./benchmark1_sigmoid 0.1 35 4 6 1 && \
    ./benchmark1_tanh 0.1 25 4 6 1 && \

    ./benchmark2_relu 0.2 10 4 6 1 && \
    ./benchmark2_relu_tanh 0.2 10 4 6 1 && \
    ./benchmark2_sigmoid 0.2 10 4 6 1 && \
    ./benchmark2_tanh 0.2 10 4 6 1 && \

    ./benchmark3_relu 0.1 60 4 6 1 && \
    ./benchmark3_relu_sigmoid 0.1 60 4 6 1 && \
    ./benchmark3_sigmoid 0.1 60 4 6 1 && \
    ./benchmark3_tanh 0.1 60 4 6 1 && \

    ./benchmark4_relu 0.02 10 4 6 1 && \
    ./benchmark4_relu_tanh 0.02 10 4 6 1 && \
    ./benchmark4_sigmoid 0.02 10 4 6 1 && \
    ./benchmark4_tanh 0.02 10 4 6 1 && \

    ./benchmark5_relu 0.02 10 4 6 1 && \
    ./benchmark5_relu_tanh 0.02 10 4 6 1 && \
    ./benchmark5_sigmoid 0.02 10 4 6 1 && \
    ./benchmark5_tanh 0.02 10 4 6 1 && \

    ./benchmark6_relu 0.02 10 4 6 1 && \
    ./benchmark6_relu_tanh 0.02 10 4 6 1 && \
    ./benchmark6_sigmoid 0.02 10 4 6 1 && \
    ./benchmark6_tanh 0.02 10 4 6 1
#    ./nn_attitude_control_sigmoid
