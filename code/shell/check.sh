for file in sim.*.tre; do
  i=$(echo "$file" | grep -o '[0-9]\+') # Extract the number from the filename
  if [ ! -f "sim.${i}.dat.csv" ]; then
    echo "Missing: sim.${i}.dat.csv"
  fi
  if [ ! -f "sim.${i}.extant.dat.csv" ]; then
    echo "Missing: sim.${i}.extant.dat.csv"
  fi
  if [ ! -f "sim.${i}.labels.csv" ]; then
    echo "Missing: sim.${i}.labels.csv"
  fi
done