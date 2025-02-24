import sys
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import random as sparse_random
from scipy.special import erfc
from copy import deepcopy

def read_matrix(filename):
    return pd.read_csv(filename, header=None).values 

def encode_ldpc(G, msg):
    return (G @ msg) % 2

def binary_symmetric_noise(codeword, flip_porb):
    noise = np.random.rand(*codeword.shape) < flip_porb
    return (codeword + noise) % 2

def decode_ldpc(H, received_bits, max_iter=100):
    for _ in range(max_iter):
        syndrome = (H @ received_bits) % 2
        if not np.any(syndrome):
            break
        flip_indices = np.argmax(H.T @ syndrome, axis=0)
        received_bits[flip_indices] ^= 1 
    return received_bits

# BER & FER Simulation
def simulate_ldpc(G, H, num_frames=1000, EbN0_range=np.arange(0, 6, 1)):
    n, k = G.shape
    ber_results = []
    fer_results = []

    for EbN0_dB in EbN0_range:
        total_bit_errors = 0
        total_frame_errors = 0
        total_bits = 0

        for _ in range(num_frames):
            msg = np.random.randint(0, 2, k)  # Random message
            print(msg)
            codeword = encode_ldpc(G, msg)
            errorenous_codeword = binary_symmetric_noise(codeword, 0.1)
            decoded = decode_ldpc(H, errorenous_codeword)

            # Count bit and frame errors
            bit_errors = np.sum(decoded != codeword)
            frame_error = bit_errors > 0

            total_bit_errors += bit_errors
            total_frame_errors += frame_error
            total_bits += len(codeword)

        # Compute BER & FER
        ber_results.append(total_bit_errors / total_bits)
        fer_results.append(total_frame_errors / num_frames)

        print(f"SNR {EbN0_dB} dB -> BER: {ber_results[-1]:.5f}, FER: {fer_results[-1]:.5f}")

    return EbN0_range, ber_results, fer_results

def depolarizing_noise(codeword, p):

    flip_prob = (2/3) *p
    noise = np.random.rand(*codeword.shape) < flip_prob
    return np.mod(codeword + noise, 2) 

# BER & FER Simulation with Depolarizing Noise
def simulate_ldpc(G, H, num_frames=1000, p_range = np.arange(0, 0.11, 0.001)):
    n, k = G.shape
    ber_results = []
    fer_results = []

    for p in p_range:
        total_bit_errors = 0
        total_frame_errors = 0
        total_bits = 0


        for _ in range(num_frames):
            msg = np.random.randint(0, 2, k)  # Random message
            codeword = encode_ldpc(G, msg)
            errorenous_codeword = depolarizing_noise(codeword, p)
            decoded = decode_ldpc(H, errorenous_codeword)

            # Count bit and frame errors
            bit_errors = np.sum(decoded != codeword)
            frame_error = bit_errors > 0

            total_bit_errors += bit_errors
            total_frame_errors += frame_error
            total_bits += len(codeword)

        # Compute BER & FER
        ber_results.append(total_bit_errors / total_bits)
        fer_results.append(total_frame_errors / num_frames)

    return p_range, ber_results, fer_results

def extract_code_config(filename):
    match = re.findall(r'(\d+)_(\d+)', filename)
    if match:
        return int(match[0][0]), int(match[0][1])


scenarios = (len(sys.argv) - 1) // 2
p_range = np.arange(0, 0.11, 0.001)
ber_results = []
fer_results = []
code_config = []
for i in range(scenarios):
    code_config.append(extract_code_config(sys.argv[i*2 + 1]))
    G = read_matrix(sys.argv[i*2 + 1])
    H = read_matrix(sys.argv[i*2 + 2])
    p_range, ber_result, fer_result = simulate_ldpc(G, H, 1000, p_range)

    ber_results.append(ber_result)
    fer_results.append(fer_result)

plt.figure(figsize=(10, 6))
for i in range(len(ber_results)): 
    n_k = code_config[i][0]
    n = code_config[i][1]
    k = n - n_k
    plt.plot(p_range, ber_results[i], label=f"BER [{n}, {k}]")

plt.xlabel('Error Probability')
plt.ylabel('Bit Error Rate (BER)')
plt.yscale('log')
plt.title('BER vs Error Probability')
plt.grid(True)
plt.legend()
plt.show()


plt.figure(figsize=(10, 6))
for i in range(len(fer_results)):
    n_k = code_config[i][0]
    n = code_config[i][1]
    k = n - n_k
    plt.plot(p_range, fer_results[i], label=f"FER [{n}, {k}]")
plt.xlabel('Error Probability')
plt.ylabel('Frame Error Rate (FER)')
plt.title('FER vs Error Probability')
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.show()