#!/usr/bin/env python3
"""Hotlink TCP proxy smoke test.

Connects to TCP proxy ports on both sides of a channel and verifies
bidirectional data transfer through the hotlink WebSocket relay.
"""

import os
import socket
import sys
import threading
import time

PORT_C1 = int(os.environ.get("PORT_CLIENT1", "9001"))
PORT_C2 = int(os.environ.get("PORT_CLIENT2", "10001"))
PAYLOAD_SIZE = int(os.environ.get("PAYLOAD_SIZE", "10000"))
TIMEOUT = int(os.environ.get("TCP_TIMEOUT", "15"))

failed = False


def recv_exact(sock, nbytes, timeout=TIMEOUT):
    sock.settimeout(timeout)
    data = b""
    while len(data) < nbytes:
        chunk = sock.recv(min(nbytes - len(data), 65536))
        if not chunk:
            raise ConnectionError(
                f"Connection closed after {len(data)}/{nbytes} bytes"
            )
        data += chunk
    return data


def connect_with_retry(host, port, retries=30, delay=0.5):
    for i in range(retries):
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.settimeout(5)
            s.connect((host, port))
            return s
        except (ConnectionRefusedError, OSError):
            if i < retries - 1:
                time.sleep(delay)
            else:
                raise
    raise RuntimeError(f"Could not connect to {host}:{port}")


def fail(msg):
    global failed
    failed = True
    print(f"  FAIL: {msg}", file=sys.stderr)


# --- Connect to both TCP proxy ports ---
print(f"Connecting to TCP proxy ports: c1={PORT_C1} c2={PORT_C2}")

sock1 = connect_with_retry("127.0.0.1", PORT_C1)
print(f"  Connected to client1 proxy port {PORT_C1}")

sock2 = connect_with_retry("127.0.0.1", PORT_C2)
print(f"  Connected to client2 proxy port {PORT_C2}")

# Let both daemons map the TCP writers
time.sleep(2)

# --- Test 1: client1 -> client2 (unidirectional) ---
print(f"\nTest 1: client1 -> client2 ({PAYLOAD_SIZE} bytes)")
payload_a = b"A" * PAYLOAD_SIZE
result_1 = [None]
err_1 = [None]


def recv_on_c2():
    try:
        result_1[0] = recv_exact(sock2, PAYLOAD_SIZE)
    except Exception as e:
        err_1[0] = str(e)


t = threading.Thread(target=recv_on_c2)
t.start()
sock1.sendall(payload_a)
t.join(timeout=TIMEOUT + 2)

if err_1[0]:
    fail(f"recv error: {err_1[0]}")
elif result_1[0] == payload_a:
    print(f"  PASS: received {len(result_1[0])} bytes correctly")
else:
    got = len(result_1[0]) if result_1[0] else 0
    fail(f"data mismatch: expected {PAYLOAD_SIZE} bytes, got {got}")

# --- Test 2: client2 -> client1 (unidirectional) ---
print(f"\nTest 2: client2 -> client1 ({PAYLOAD_SIZE} bytes)")
payload_b = b"B" * PAYLOAD_SIZE
result_2 = [None]
err_2 = [None]


def recv_on_c1():
    try:
        result_2[0] = recv_exact(sock1, PAYLOAD_SIZE)
    except Exception as e:
        err_2[0] = str(e)


t = threading.Thread(target=recv_on_c1)
t.start()
sock2.sendall(payload_b)
t.join(timeout=TIMEOUT + 2)

if err_2[0]:
    fail(f"recv error: {err_2[0]}")
elif result_2[0] == payload_b:
    print(f"  PASS: received {len(result_2[0])} bytes correctly")
else:
    got = len(result_2[0]) if result_2[0] else 0
    fail(f"data mismatch: expected {PAYLOAD_SIZE} bytes, got {got}")

# --- Test 3: bidirectional simultaneous ---
print(f"\nTest 3: bidirectional ({PAYLOAD_SIZE} bytes each way)")
payload_c = b"C" * PAYLOAD_SIZE
payload_d = b"D" * PAYLOAD_SIZE
result_3c1 = [None]
result_3c2 = [None]
err_3 = []


def bidir_c1():
    try:
        sock1.sendall(payload_c)
        result_3c1[0] = recv_exact(sock1, PAYLOAD_SIZE)
    except Exception as e:
        err_3.append(f"c1: {e}")


def bidir_c2():
    try:
        sock2.sendall(payload_d)
        result_3c2[0] = recv_exact(sock2, PAYLOAD_SIZE)
    except Exception as e:
        err_3.append(f"c2: {e}")


t1 = threading.Thread(target=bidir_c1)
t2 = threading.Thread(target=bidir_c2)
t1.start()
t2.start()
t1.join(timeout=TIMEOUT + 2)
t2.join(timeout=TIMEOUT + 2)

if err_3:
    fail(f"bidir errors: {'; '.join(err_3)}")
elif result_3c2[0] == payload_c and result_3c1[0] == payload_d:
    print(f"  PASS: bidirectional transfer correct")
else:
    details = []
    if result_3c2[0] != payload_c:
        got = len(result_3c2[0]) if result_3c2[0] else 0
        details.append(f"c2 got {got}/{PAYLOAD_SIZE}")
    if result_3c1[0] != payload_d:
        got = len(result_3c1[0]) if result_3c1[0] else 0
        details.append(f"c1 got {got}/{PAYLOAD_SIZE}")
    fail(f"bidir mismatch: {'; '.join(details)}")

# --- Test 4: larger payload ---
BIG = PAYLOAD_SIZE * 10
print(f"\nTest 4: large payload ({BIG} bytes each way)")
payload_e = bytes(range(256)) * (BIG // 256 + 1)
payload_e = payload_e[:BIG]
payload_f = bytes(reversed(range(256))) * (BIG // 256 + 1)
payload_f = payload_f[:BIG]
result_4c1 = [None]
result_4c2 = [None]
err_4 = []


def big_c1():
    try:
        sock1.sendall(payload_e)
        result_4c1[0] = recv_exact(sock1, BIG)
    except Exception as e:
        err_4.append(f"c1: {e}")


def big_c2():
    try:
        sock2.sendall(payload_f)
        result_4c2[0] = recv_exact(sock2, BIG)
    except Exception as e:
        err_4.append(f"c2: {e}")


t1 = threading.Thread(target=big_c1)
t2 = threading.Thread(target=big_c2)
t1.start()
t2.start()
t1.join(timeout=TIMEOUT * 2)
t2.join(timeout=TIMEOUT * 2)

if err_4:
    fail(f"large payload errors: {'; '.join(err_4)}")
elif result_4c2[0] == payload_e and result_4c1[0] == payload_f:
    print(f"  PASS: large bidirectional transfer correct")
else:
    details = []
    if result_4c2[0] != payload_e:
        got = len(result_4c2[0]) if result_4c2[0] else 0
        details.append(f"c2 got {got}/{BIG}")
    if result_4c1[0] != payload_f:
        got = len(result_4c1[0]) if result_4c1[0] else 0
        details.append(f"c1 got {got}/{BIG}")
    fail(f"large bidir mismatch: {'; '.join(details)}")

# --- Cleanup ---
sock1.close()
sock2.close()

if failed:
    print("\nSOME TESTS FAILED", file=sys.stderr)
    sys.exit(1)
else:
    print("\nAll hotlink TCP proxy tests passed!")
    sys.exit(0)
