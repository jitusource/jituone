import sys, string
import os, commands, re, glob
import csv

root_total=10466

print 'ROOT probes total:', root_total

F3= csv.DictReader(open("ROOT_match-3.5tr_pc0.98_fracL0.5.csv"))#10466 
shoot_unq3=[]
f3D=[]
for seq in F3: #[{'#Probesets': 'Mtr.2632.1.S1_at', 'Representative Public ID': 'BI311277'}]
    probe=seq['Probe'].strip()
    shoot_unq3.append(probe)
    f3D.append(seq)

print '# found in v3.5: ', len(shoot_unq3)

shoot_unq3=list(set(shoot_unq3))

F4= csv.DictReader(open("ROOT_match-4_pct0.98_fracL0.5.csv"))#10466 
shoot_unq4=[]
f4D=[]
for seq in F4: #[{'#Probesets': 'Mtr.2632.1.S1_at', 'Representative Public ID': 'BI311277'}]
    probe=seq['Probe'].strip()
    shoot_unq4.append(probe)
    f4D.append(seq)
     

print '# found in v4.0: ', len(shoot_unq3)
shoot_unq4=list(set(shoot_unq4))###




print 'uniq found in v3.5: ', len(shoot_unq3)
print 'uniq found in v4.0: ', len(shoot_unq4)

print 'in model 4  but not in 3.5v:', len( list(set(shoot_unq4) - set(shoot_unq3)) )
print 'in 3.5v3, but not in 4:', len( list(set(shoot_unq3) - set(shoot_unq4)) )
print 'Common to boths: ', len(set(shoot_unq3) & set(shoot_unq4) )# common in both lists

master=list(set(shoot_unq3 + shoot_unq4))
print 'Merging both:', len(master)
print 'Coverage percentage after merging:', 100*len(master)/float(root_total)


out=open('ROOT_score3n4_' + '_pct0.98_fracL0.5.csv','w')

COL=['Probe', 'Public_ID', 'Description', 'Title', 'where', 'score' ] 
top_tile=','.join(COL) + '\n'
out.write(top_tile)

for mas in master: #[{'#Probesets': 'Mtr.2632.1.S1_at', 'Representative Public ID': 'BI311277'}]
    max_score=-1
    max_seq=None
    for seq in f3D:
    	probe=seq['Probe'].strip()
    	if mas == probe:
       	   score=float(seq['score'].strip())
           #print score, seq
           if score > max_score:
           	   max_score=score
           	   max_seq=seq
           	   max_seq['where']='v3.5'

    for seq in f4D:
    	probe=seq['Probe'].strip()
    	if mas == probe:
       	   score=float(seq['score'].strip())
           #print score, seq
           if score > max_score:
           	   max_score=score
           	   max_seq=seq
           	   max_seq['where']='v4.0'

    #print '----max ---'
    if max_score  > -1:
       #print max_score
       #print max_seq
       top_tile=','.join([str(max_seq[k]) for k in COL]) + '\n' 
       out.write(top_tile)

    #if mas == 
    #best_hits= { 'Probe': probe, 'Public_ID': pub_id, 'Title': full_title, 'Description': descript, 'Match_serial' : str(dups), 'score':  str(bit_score) } ) ###'Full_length_Duplicate_match': str(dups),  'Hit_100_pctID':str(hit_hundred)})
                           

print ' DONE '

