for f in `cat run.host`; do echo $f; ./measure-binary.sh $f; done
a=`grep flags /proc/cpuinfo | tail -1`
if [[$a == *"avx2"*]] then
  ./measure-binary.sh lmp_intel_eval_Cvi;
  ./measure-binary.sh lmp_intel_eval_Cvj;
elif [[$a == *"avx"*]] then
  ./measure-binary.sh lmp_intel_eval_Avi;
  ./measure-binary.sh lmp_intel_eval_Avj;
elif [[$a == *"sse4_2"*]] then
  ./measure-binary.sh lmp_intel_eval_Svi;
  ./measure-binary.sh lmp_intel_eval_Svj;
fi;
/usr/bin/ssh `hostname -s`-mic0 "cd ${PWD##$HOME/}; for f in \`cat run.mic\`; do echo \$f; ./measure-binary.sh \$f; done"

