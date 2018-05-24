no version (v1)
target: signif peaks in ensg_id and at least one motif in signif peaks
KD gene list in K562 or Hepg2 (in one of them is enough)
This results in 151 target genes 
non target: same len as target
151, 151

v2:
target: same as v1
non target: larger non target list  (if 3 times is in effect)
153, 423

v3:
target: same as v1
non target: larger than non target list (removed if 3 times but selected targets from bins that don’t have any non target)
160
295


v4: 
target: signify peak, motif in signify peak, not KD on the other side
non target: one to one match
1046-1046



v5:
target: signif peak + motif in signif peak, KD threshold decreased to LFC:0 meaning fold change is 1
non target: don’t sample from non target list when non targets are larger than targets, choose non targets that are longest…
multfactor: 1.2 

T: 423
N: 440
binsize: 120

v6: binsize: 100
targets: same as v5:
non target: same as v5: multfactor:1.1
T: 433
N: 413

PTBP1
default code couldn't create length matched background set, I set select_oneortwo_targets to False
check the code
p-value mannwhitneyu in terms of length  PTBP1 0.40373825913 1227.0 1176.5
PTBP1 215 226

RBFOX2:
KD check targets removed it
length matched problem, select_oneortwo_targets = False


IGF2BP2
num_first_bins  13 8
p-value mannwhitneyu in terms of length  IGF2BP2 0.0583052651635 1583.0 1433.0
IGF2BP2 210 201

the numbers (T original is of size ~400) could go up if non target list gets larger.


TAF15
removed motif_check and KD_check from target completely
138 133

KHSRP:
got rid of KD check for both target and non target
select_oneortwo_targets = False
KHSRP 413 444

HNRNPK
fine
290 313

U2AF2
select_oneortwo_targets = False
U2AF2 141 149

HNRNPC
removed KD check from non target
select_oneortwo_targets = False
HNRNPC 122 127





RBM22 
removed KD check from targets
removed motif signify peaks from targets
RBM22 267 281



TARDBP 
select_oneortwo_targets = False
TARDBP 300 309


QKI 
select_oneortwo_targets = False
QKI 206 223


SRSF1 
fine
SRSF1 256 255




PCBP2 
select_oneortwo_targets = False
PCBP2 687 762


HNRNPM 
removed KD check, removed motif in all peaks as well, motif check only counts # of motifs
HNRNPM 135 123


EIF4G2 
removed KD check, removed motif in all peaks as well, motif check only counts # of motifs
select_oneortwo_targets = False
EIF4G2 270 292


TIA1
select_oneortwo_targets = False
TIA1 132 134



FMR1
select_oneortwo_targets = False
320 348

FXR1 
removed motif check and KD completely from target
FXR1 479 514

IGF2BP1 
select_oneortwo_targets = False
IGF2BP1 177 185



FXR2 
removed KD check from target
FXR2 133 140



EWSR1
removed KD check and signify peak motif from target
EWSR1 105 107


SRSF9
removed KD and motif check from non target
removed KD check from target
SRSF9 221 243


TRA2A 
removed KD check from target
TRA2A 164 172



SFPQ 2300 2300
select_oneortwo_targets = False
SFPQ 64 65
removed motif check and KD check completely from target list

SRSF7
removed KD and motif check from target list
SRSF7 514 579


FUBP3 
select_one False

FUBP3 150 152


KHDRBS1 
KD check and motif in peak check is gone from target list
motif_check > 0 remains only
KD check is removed from non target list
KHDRBS1 112 99

IGF2BP3 
select_one = False
IGF2BP3 154 162


HNRNPU 
removed motif check and KD check from target
KD check turned to or for non targets
HNRNPU 126 132



HNRNPA1 
select_one = False
KD check removed from target
motif check all peaks remains
KD check turned to or for non targets

HNRNPA1 214 230


PTBP1 7133 7133 why two of these?


HNRNPF 
select_ False
motif_check >0 remains for target
removed KD and motifs in peaks from target
HNRNPF 130 138


LIN28A 
removed KD check from target
LIN28A 137 133


FUS 
select_one = False
179 188


HNRNPD 
removed KD check from targets
select = False
HNRNPD 381 422


RBM47 
no p-value on piranha based CLIP 
so took all 
changed KD to or
RBM47 101 105


HNRNPA2B1
removed all KD check and motif check from targets
select_one = False
107 111


HNRNPL 
select_one = False
score based eCLIP threshold 1
HNRNPL 221 239

RBM5 
select = False
removed KD check from target list
score based eCLIP threshold 1
RBM5 386 438

ELAVL1:
select = False
used Mukharjee_dataset as KD
711 738

1) remove KD check from non target, first make it or then remove it 
2) remove motif check from non target
3) remove KD check from target, 
4) change motif check to all_peaks 
5) change motif check to anywhere i.e., dict_motif[esng_id]
6) remove motif check completely from target





