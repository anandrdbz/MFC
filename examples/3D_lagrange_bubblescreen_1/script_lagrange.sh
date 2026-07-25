#!/bin/bash
#!/bin/bash

cd ../
for i in {2..40}; do
    rm -rf "3D_lagrange_qbmm_${i}"
    mkdir -p "3D_lagrange_qbmm_${i}"
    cd "3D_lagrange_qbmm_${i}" || exit 1

    cp ../3D_lagrange_qbmm_1/main.py ./
    cp ../3D_lagrange_qbmm_1/case.py ./
    cp ../3D_lagrange_qbmm_1/py_functions.py ./

    mkdir -p input

    sed -i "s/rnd\.seed(1)/rnd.seed(${i})/g" py_functions.py

    python3 main.py

    cd ../../ || exit 1

    . ./mfc.sh load -c f -m g

    ./mfc.sh run "examples/3D_lagrange_qbmm_${i}/case.py" \
        -j 8 -N 1 -n 8 \
        --case-optimization \
        -c frontier \
        -t pre_process simulation \
        --no-gpu \
        -e batch \
        -w 01:00:00 \
        -a cfd154

    cd examples || exit 1
done
