for i in 0 1 2 3 4 5 6 7 8 9; do
  for j in si-bulk cnt; do
    ./bin/$1 -in in.$j -pk intel 0 mode double -sf intel > out/`hostname -s`.${j}.double.${1}.$i
    ./bin/$1 -in in.$j -pk intel 0 mode single -sf intel > out/`hostname -s`.${j}.single.${1}.$i
    ./bin/$1 -in in.$j -pk intel 0 mode mixed -sf intel > out/`hostname -s`.${j}.mixed.${1}.$i
    [[ $1 =~ lmp_intel_eval_.vi ]] && ./bin/$1 -in in.$j > out/`hostname -s`.${j}.original.${1}.$i
  done
done
