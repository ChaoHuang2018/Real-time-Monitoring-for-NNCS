make && \
    # ./benchmark1_relu 0.01 35 4 6 && \  # empty poly
    ./benchmark1_relu_tanh 0.01 35 4 6 && \
    ./benchmark1_sigmoid 0.01 35 4 6 && \
    ./benchmark1_tanh 0.01 35 4 6 && \
    ./benchmark2_relu 0.01 10 4 6 && \
    ./benchmark2_relu_tanh 0.01 10 4 6 && \
    ./benchmark2_sigmoid 0.01 10 4 6 && \
    # ./benchmark2_tanh 0.01 10 4 6 && \  # dimension mismatch
    ./benchmark3_relu 0.01 60 4 6 && \
    ./benchmark3_relu_sigmoid 0.01 60 4 6 && \
    ./benchmark3_sigmoid 0.01 60 4 6 && \
    ./benchmark3_tanh 0.01 60 4 6 && \
    ./benchmark4_relu 0.01 10 4 6 && \
    ./benchmark4_relu_tanh 0.01 10 4 6 && \
    ./benchmark4_sigmoid 0.01 10 4 6 && \
    ./benchmark4_tanh 0.01 10 4 6 && \
    # ./benchmark5_relu 0.01 10 4 6 && \  # empty poly
    ./benchmark5_relu_tanh 0.01 10 4 6 && \
    ./benchmark5_sigmoid 0.01 10 4 6 && \
    ./benchmark5_tanh 0.01 10 4 6 && \
    ./benchmark6_relu 0.01 10 4 6 && \
    # ./benchmark6_relu_tanh 0.01 10 4 6 && \  # empty poly
    ./benchmark6_sigmoid 0.01 10 4 6 && \
    ./benchmark6_tanh 0.01 10 4 6
#    ./nn_attitude_control_sigmoid
