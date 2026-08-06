#!/bin/bash
for np in 1 2 4 6 8 10 12; do
    echo "Aguardando 30 segundos para esfriar (np=$np)..."
    sleep 30
    echo "=== Testando np=$np ==="
    MSCD_GPU=1 ./baseline/regressao-gpu.sh $np
    echo ""
done
