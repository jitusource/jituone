import sys, string
import os, commands, re, glob
import csv
from Bio.Blast import NCBIStandalone
blast_parser = NCBIStandalone.BlastParser()
  

def sanitize(name):
    acceptable_characters = string.ascii_letters + string.digits + '_' + ' ' + '|''|'
    new_name = ''.join(character for character in name if character in acceptable_characters)
    assert len(new_name) > 0 and new_name[0] not in string.digits+'_', "change %r to %r; is still invalid" % (name, new_name)
    return new_name



### step 1 ### get all the names probe Id and pub_Id names from CSV and remove if duplicate names in Probes #
dRoot= csv.DictReader(open("shoots_probe_set.csv"))#10466

uniq_probe_set=[]
many=0
PROBE_PUB=[]
for seq in dRoot:  
    probe=seq['#Probesets'].strip()
    pub_id=seq['Representative Public ID'].strip()

    if probe  in uniq_probe_set:
       print 'duplicate probe in CSV so bumping it : ', probe
       continue
    else:
       uniq_probe_set.append(probe) 
       PROBE_PUB.append([probe,pub_id]) 

print len(PROBE_PUB)



#for fraction_match in [0.35, 0.55, 0.65, 0.75, 0.85,  0.95, 0.999]:
for fraction_match in [0.50]:
    pct_id_thresh=0.98
    #fraction_match=0.75
    out=open('SHOOT_match_3.5tr_pct' + str(pct_id_thresh) +'_fracL'+ str(fraction_match )  + '.csv','w')
    #out=open('full_matches_roots.csv','w')

    COL=['Probe', 'Public_ID', 'Description', 'Title','Match_serial', 'score' ] 
    top_tile=','.join(COL) + '\n'
    out.write(top_tile)

    uniq_probe_found=[]
    ### step 2 ### Collect  number of hits in list 
    for probe, pub_id in PROBE_PUB:
        try:
         result_handle = open(probe + '.fasta.blast')
        except:
          continue 
        blast_iterator = NCBIStandalone.Iterator(result_handle, blast_parser)       
        hits=[]
        dups=0
        hit_hundred=0
        for blast_record in blast_iterator:
            query_len=blast_record.query_letters ## how many letters in query     
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    #print hsp.score  # bracket 
                    bit_score=hsp.bits #993
                    #print probe + '.fasta.blast'
                    perc=hsp.identities
                    if (float(perc[0])/float(perc[1])) > pct_id_thresh :### hundred percent identity
                       full_len_match=False 
                       align_title=sanitize(alignment.title.strip())
                       align_title.replace(",", "")# remove commas 
                       clean_title=' '.join(align_title.split()) 
                       full_title=clean_title[:]
                       descript=clean_title.split('|')[1]
                       hit_hundred+=1                   
                       #----------------------------------------------------
                       if (float(perc[1])/float(query_len)) > fraction_match:### if its full length match or not  ###
                          full_len_match=True
                          dups+=1 
                          hits.append({ 'Probe': probe, 'Public_ID': pub_id, 'Title': full_title, 'Description': descript, 'Match_serial' : str(dups), 'score':  str(bit_score) } ) ###'Full_length_Duplicate_match': str(dups),  'Hit_100_pctID':str(hit_hundred)})
                          if probe not in uniq_probe_found:
                             uniq_probe_found.append(probe)

        ### write hits with only full length matches #####
        #COL=['Probe', 'Public_ID', 'Description', 'Title', 'Full_length_Duplicate_match', 'Hit_100_pctID'] 
        for topD in hits:
            top_tile=','.join([str(topD[k]) for k in COL]) + '\n' 
            out.write(top_tile)

    print 'Finding match coverage shoots 3.5V '
    print '# uniq_present: ',len(uniq_probe_set)
    print '# uniq_Found: ',len(uniq_probe_found)
    print '# match: percent ',  len(uniq_probe_found)/(1.*len(uniq_probe_set))
    print 'percent ID >',  pct_id_thresh*100
    print 'Fraction Match >',  fraction_match*100
    print  'ROOTS DONE '
    print '------------------------------------'


print 'DONE: this FULL path =', os.path.abspath(os.path.dirname(sys.argv[0]))  + '/'+  sys.argv[0]
print '---- done '
