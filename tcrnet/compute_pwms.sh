MIRSCRIPT="java -cp mir-1.0-SNAPSHOT.jar com.milaboratory.mir.scripts.Examples"

$MIRSCRIPT clonotype-summary-stats -S Human -G TRA -U AA_CDR3_PWM -F VDJtools -I human.tra.aa.txt -O pwms/human.tra.aa
$MIRSCRIPT clonotype-summary-stats -S Human -G TRB -U AA_CDR3_PWM -F VDJtools -I human.trb.aa.txt -O pwms/human.trb.aa
$MIRSCRIPT clonotype-summary-stats -S Mouse -G TRA -U AA_CDR3_PWM -F VDJtools -I mouse.tra.aa.txt -O pwms/mouse.tra.aa
$MIRSCRIPT clonotype-summary-stats -S Mouse -G TRB -U AA_CDR3_PWM -F VDJtools -I mouse.trb.aa.txt -O pwms/mouse.trb.aa