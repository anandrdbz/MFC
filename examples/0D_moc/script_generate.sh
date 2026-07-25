#Generate
for i in $(seq 1 1000); do
    CASE_DIR="0D_moc_${i}"
    mkdir -p "$CASE_DIR"
    cp case.py "${CASE_DIR}/case.py"

    sed -i '' \
        "s/^sample_index[[:space:]]*=[[:space:]]*[0-9][0-9]*/sample_index = ${i}/" \
        "${CASE_DIR}/case.py"
done
