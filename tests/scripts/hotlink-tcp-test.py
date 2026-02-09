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
BURST_COUNT = int(os.environ.get("BURST_COUNT", "0"))
BURST_ONLY = os.environ.get("BURST_ONLY", "0") == "1"
SOCK_BUF = int(os.environ.get("SOCK_BUF_BYTES", "0"))

failed = False
results = []


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
            if SOCK_BUF > 0:
                s.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, SOCK_BUF)
                s.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, SOCK_BUF)
            # Use blocking mode for sustained throughput tests.
            s.settimeout(None)
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


def record_result(name, bytes_total, seconds, ok):
    mbps = (bytes_total / (1024 * 1024)) / seconds if seconds > 0 else 0.0
    results.append(
        {
            "name": name,
            "bytes": int(bytes_total),
            "seconds": float(seconds),
            "mbps": float(mbps),
            "ok": bool(ok),
        }
    )


def run_burst_unidirectional(sock_tx, sock_rx, payload_size, count, timeout):
    """Send N framed payloads as fast as possible in one direction."""
    payload = b"Z" * payload_size
    header = payload_size.to_bytes(4, "big")
    frame = header + payload
    got = [0]
    err = [None]
    start = time.perf_counter()

    def recv_all():
        try:
            for _ in range(count):
                size = int.from_bytes(recv_exact(sock_rx, 4, timeout=timeout), "big")
                if size != payload_size:
                    raise ValueError(f"frame size mismatch: expected {payload_size}, got {size}")
                recv_exact(sock_rx, size, timeout=timeout)
                got[0] += 1
        except Exception as e:
            err[0] = str(e)

    t = threading.Thread(target=recv_all)
    t.start()
    for _ in range(count):
        sock_tx.sendall(frame)
    t.join(timeout=max(timeout, 10) + count)
    elapsed = time.perf_counter() - start
    ok = err[0] is None and got[0] == count
    return ok, got[0], elapsed, err[0]


# --- Connect to both TCP proxy ports ---
print(f"Connecting to TCP proxy ports: c1={PORT_C1} c2={PORT_C2}")

sock1 = connect_with_retry("127.0.0.1", PORT_C1)
print(f"  Connected to client1 proxy port {PORT_C1}")

sock2 = connect_with_retry("127.0.0.1", PORT_C2)
print(f"  Connected to client2 proxy port {PORT_C2}")

# Let both daemons map the TCP writers
time.sleep(2)

if BURST_ONLY:
    if BURST_COUNT <= 0:
        fail("BURST_ONLY=1 requires BURST_COUNT > 0")
        print("\nSOME TESTS FAILED", file=sys.stderr)
        sys.exit(1)
    print(
        f"\nBurst benchmark: client1 -> client2, "
        f"{BURST_COUNT} x {PAYLOAD_SIZE} bytes"
    )
    ok, got_count, elapsed, err = run_burst_unidirectional(
        sock1, sock2, PAYLOAD_SIZE, BURST_COUNT, TIMEOUT
    )
    total_bytes = got_count * PAYLOAD_SIZE
    avg_ms = (elapsed / got_count * 1000.0) if got_count > 0 else 0.0
    throughput = (total_bytes / (1024 * 1024)) / elapsed if elapsed > 0 else 0.0
    status = "PASS" if ok else "FAIL"
    print(
        f"  burst-c1->c2  {status}  files={got_count}/{BURST_COUNT}  "
        f"file_size={PAYLOAD_SIZE}  total_bytes={total_bytes}  "
        f"time={elapsed:.4f}s  avg_file={avg_ms:.3f}ms  throughput={throughput:.2f} MiB/s"
    )
    if err:
        fail(f"burst error: {err}")
    sock1.close()
    sock2.close()
    if failed:
        print("\nSOME TESTS FAILED", file=sys.stderr)
        sys.exit(1)
    print("\nBurst benchmark completed successfully")
    sys.exit(0)

# --- Test 1: client1 -> client2 (unidirectional) ---
print(f"\nTest 1: client1 -> client2 ({PAYLOAD_SIZE} bytes)")
payload_a = b"A" * PAYLOAD_SIZE
result_1 = [None]
err_1 = [None]
test1_start = time.perf_counter()


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
    record_result("c1->c2", PAYLOAD_SIZE, time.perf_counter() - test1_start, False)
elif result_1[0] == payload_a:
    print(f"  PASS: received {len(result_1[0])} bytes correctly")
    record_result("c1->c2", PAYLOAD_SIZE, time.perf_counter() - test1_start, True)
else:
    got = len(result_1[0]) if result_1[0] else 0
    fail(f"data mismatch: expected {PAYLOAD_SIZE} bytes, got {got}")
    record_result("c1->c2", got, time.perf_counter() - test1_start, False)

# --- Test 2: client2 -> client1 (unidirectional) ---
print(f"\nTest 2: client2 -> client1 ({PAYLOAD_SIZE} bytes)")
payload_b = b"B" * PAYLOAD_SIZE
result_2 = [None]
err_2 = [None]
test2_start = time.perf_counter()


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
    record_result("c2->c1", PAYLOAD_SIZE, time.perf_counter() - test2_start, False)
elif result_2[0] == payload_b:
    print(f"  PASS: received {len(result_2[0])} bytes correctly")
    record_result("c2->c1", PAYLOAD_SIZE, time.perf_counter() - test2_start, True)
else:
    got = len(result_2[0]) if result_2[0] else 0
    fail(f"data mismatch: expected {PAYLOAD_SIZE} bytes, got {got}")
    record_result("c2->c1", got, time.perf_counter() - test2_start, False)

# --- Test 3: bidirectional simultaneous ---
print(f"\nTest 3: bidirectional ({PAYLOAD_SIZE} bytes each way)")
payload_c = b"C" * PAYLOAD_SIZE
payload_d = b"D" * PAYLOAD_SIZE
result_3c1 = [None]
result_3c2 = [None]
err_3 = []
test3_start = time.perf_counter()


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
    record_result("bidir-small", PAYLOAD_SIZE * 2, time.perf_counter() - test3_start, False)
elif result_3c2[0] == payload_c and result_3c1[0] == payload_d:
    print(f"  PASS: bidirectional transfer correct")
    record_result("bidir-small", PAYLOAD_SIZE * 2, time.perf_counter() - test3_start, True)
else:
    details = []
    if result_3c2[0] != payload_c:
        got = len(result_3c2[0]) if result_3c2[0] else 0
        details.append(f"c2 got {got}/{PAYLOAD_SIZE}")
    if result_3c1[0] != payload_d:
        got = len(result_3c1[0]) if result_3c1[0] else 0
        details.append(f"c1 got {got}/{PAYLOAD_SIZE}")
    fail(f"bidir mismatch: {'; '.join(details)}")
    got = (len(result_3c1[0]) if result_3c1[0] else 0) + (len(result_3c2[0]) if result_3c2[0] else 0)
    record_result("bidir-small", got, time.perf_counter() - test3_start, False)

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
test4_start = time.perf_counter()


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
    record_result("bidir-large", BIG * 2, time.perf_counter() - test4_start, False)
elif result_4c2[0] == payload_e and result_4c1[0] == payload_f:
    print(f"  PASS: large bidirectional transfer correct")
    record_result("bidir-large", BIG * 2, time.perf_counter() - test4_start, True)
else:
    details = []
    if result_4c2[0] != payload_e:
        got = len(result_4c2[0]) if result_4c2[0] else 0
        details.append(f"c2 got {got}/{BIG}")
    if result_4c1[0] != payload_f:
        got = len(result_4c1[0]) if result_4c1[0] else 0
        details.append(f"c1 got {got}/{BIG}")
    fail(f"large bidir mismatch: {'; '.join(details)}")
    got = (len(result_4c1[0]) if result_4c1[0] else 0) + (len(result_4c2[0]) if result_4c2[0] else 0)
    record_result("bidir-large", got, time.perf_counter() - test4_start, False)

print("\nPerformance summary:")
for row in results:
    status = "PASS" if row["ok"] else "FAIL"
    print(
        f"  {row['name']:<12} {status}  bytes={row['bytes']}  "
        f"time={row['seconds']:.4f}s  throughput={row['mbps']:.2f} MiB/s"
    )

if BURST_COUNT > 0:
    print(
        f"\nBurst benchmark: client1 -> client2, "
        f"{BURST_COUNT} x {PAYLOAD_SIZE} bytes"
    )
    ok, got_count, elapsed, err = run_burst_unidirectional(
        sock1, sock2, PAYLOAD_SIZE, BURST_COUNT, TIMEOUT
    )
    total_bytes = got_count * PAYLOAD_SIZE
    avg_ms = (elapsed / got_count * 1000.0) if got_count > 0 else 0.0
    throughput = (total_bytes / (1024 * 1024)) / elapsed if elapsed > 0 else 0.0
    status = "PASS" if ok else "FAIL"
    print(
        f"  burst-c1->c2  {status}  files={got_count}/{BURST_COUNT}  "
        f"file_size={PAYLOAD_SIZE}  total_bytes={total_bytes}  "
        f"time={elapsed:.4f}s  avg_file={avg_ms:.3f}ms  throughput={throughput:.2f} MiB/s"
    )
    if err:
        fail(f"burst error: {err}")

# --- Cleanup ---
sock1.close()
sock2.close()

if failed:
    print("\nSOME TESTS FAILED", file=sys.stderr)
    sys.exit(1)
else:
    print("\nAll hotlink TCP proxy tests passed!")
    sys.exit(0)
