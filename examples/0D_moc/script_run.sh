cd ../../

#Run
for i in $(seq 1 1000); do
    ./mfc.sh run examples/0D_moc/0D_moc_${i}/case.py -t pre_process simulation --case-optimization
done
