# adaptiveTBV for predicting zoonotic risk to human of a TBV RdRp CDS. 
main steps:
1. initial raw data with sequences and innotaions was under "main/data/" and with label file under "main/labels/".
2. a sampling was suggested for further analysis such as counting.
3. counting sequences for DCR/codonpair via "batch_dcrcp_count.py" under "main/counting/" from the sequence file under "main/data/".
4. unsupervised learning to overview the data distribution based on DCR/codonpair was available via "unsupervisedlearning.py" under the path of "main/unsupervised_clustering/".
5. supvervised predictor training was applicable with the script "trainer.py" under the path of "main/supervised_training/".
6. A trained model of "xx.pt" under the path of "main/supervised_training/" was available for predict the zoonotic risk of a TBV RdRp CDS post model training, and ten trained models were also uploaded as supplementary data9, in the submitted manuscript to the journal, Med (https://www.cell.com/Med/home).

Citation:
{"Title": "The diversity of emerging tick-borne viruses globally: from discoveries to zoonotic risk assessment",
  "Authors": "Mei-Qi Zhang1,9, Jing Li1,9, Xiao-Long Lv2,9, Jin-Jin Chen1, Shu-Zhen Han2, Yun-Bo Qiu1,3, Xiao-Hu Han4, Guang-Qian Si1, Zhen-Yu Hu1, Hong-Han Ge1,5, Xiao-Ai Zhang1, Chang Li5, Tao Jiang1, Max L. Mehlman7,8, Simon I. Hay7,8,*, Li-Qun Fang1,3,*, Wei Liu1,10,*",
"Affiliations": "1State Key Laboratory of Pathogen and Biosecurity, Academy of Military Medical Sciences, Beijing, 100071, China.
2 Second Affiliated Hospital of Inner Mongolia University for the Nationalities, Inner Mongolia General Forestry Hospital, Yakeshi, 022150, China.
3 School of Public Health, Guizhou Medical University, Guiyang, 550004, China.
4 Key Laboratory of Livestock Infectious Diseases, Shenyang Agricultural University, Shenyang, 110866, China.
5School of Public Health, Shandong First Medical University and Shandong Academy of Medical Sciences, Jinan, 250117, China.
6Chinese Academy of Medical Sciences, Changchun Veterinary Research Institute, Chinese Academy of Agricultural Sciences, Changchun, 130122, China.
7Department of Health Metrics Sciences, School of Medicine, University of Washington, Washington, USA. 
8Institute for Health Metrics and Evaluation, University of Washington, Washington, USA."
"Contact information": "9These authors contributed equally.10Lead contact, *Correspondence and requests for materials should be addressed to Prof Simon I. Hay (email: sihay@uw.edu), Prof Li-Qun Fang (email: fang_lq@163.com) or Prof Wei Liu (email: lwbime@163.com)."
