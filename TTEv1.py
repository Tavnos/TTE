from random import choices
from tkinter import (Canvas, Frame, Text, Entry, Label, Button, Radiobutton, Checkbutton, StringVar, Tk)


class Nucleotide_Data:
    dna_base = []
    adenine = {'char':'a','dna':'a','match':'t','rna':'a'}
    guanine = {'char':'g','dna':'g','match':'c','rna':'g'}
    cytosine = {'char':'c','dna':'c','match':'g','rna':'c'}
    thymine = {'char':'t','dna':'t','match':'a','rna':'u'}
    uracil = {'char':'u','dna':'t','match':'a','rna':'u'}
    nucleotide = (adenine,guanine,cytosine,thymine,uracil)
    def __init__(self):
        for i in range(len(self.nucleotide)-1):
            self.dna_base.append(self.nucleotide[i]['char'])
class Amino_Data:
    amino_acid = []
    methionine = {'config':('atg'),'name':'Methionine','letter_code':'M','tri_char':'Met','amino_class':'sulfuric','amino_polarity':'non-polar','amino_charge':'neutral','amino_hydropathy':1.9,'amino_weight':149.208}
    leucine = {'config':('tta','ttg','ctt','ctc','cta','ctg'),'name':'leucine','letter_code':'L','tri_char':'Leu','amino_class':'aliphatic','amino_polarity':'non-polar','amino_charge':'neutral','amino_hydropathy':3.8,'amino_weight':131.175}
    phenylalanine = {'config':('ttt','ttc'),'name':'phenylalanine','letter_code':'F','tri_char':'Phe','amino_class':'aromatic','amino_polarity':'non-polar','amino_charge':'neutral','amino_hydropathy':2.8,'amino_weight':165.192}
    serine = {'config':('tct','tcc','tca','tcg','agc','agt'),'name':'serine','letter_code':'S','tri_char':'Ser','amino_class':'hydroxylic','amino_polarity':'polar','amino_charge':'neutral','amino_hydropathy':1.9,'amino_weight':149.208}
    cysteine = {'config':('tgt','tgc'),'name':'cysteine','letter_code':'C','tri_char':'Cys','amino_class':'sulfuric','amino_polarity':'non-polar','amino_charge':'neutral','amino_hydropathy':2.5,'amino_weight':121.154}
    tryptophan = {'config':('tgg'),'name':'tryptophan','letter_code':'W','tri_char':'Trp','amino_class':'aromatic','amino_polarity':'non-polar','amino_charge':'neutral','amino_hydropathy':-0.9,'amino_weight':204.228}
    proline = {'config':('ccc','cca','cct','ccg'),'name':'proline','letter_code':'P','tri_char':'Pro','amino_class':'cyclic','amino_polarity':'non-polar','amino_charge':'neutral','amino_hydropathy':-1.6,'amino_weight':115.132}
    histidine = {'config':('cat','cac'),'name':'histidine','letter_code':'H','tri_char':'His','amino_class':'basic-aromatic','amino_polarity':'basic-polar','amino_charge':'neutral-positive','amino_hydropathy':-3.2,'amino_weight':155.156}
    glutamine = {'config':('caa','cag'),'name':'glutamine','letter_code':'Q','tri_char':'Gln','amino_class':'amide','amino_polarity':'polar','amino_charge':'neutral','amino_hydropathy':-3.5,'amino_weight':146.146}
    arginine = {'config':('cgg','cgt','cgc','cga','agg','aga'),'name':'arginine','letter_code':'R','tri_char':'Arg','amino_class':'basic','amino_polarity':'basic-polar','amino_charge':'positive','amino_hydropathy':-4.5,'amino_weight':174.203}
    isoleucine = {'config':('att','atc','ata'),'name':'isoleucine','letter_code':'I','tri_char':'Ile','amino_class':'aliphatic','amino_polarity':'non-polar','amino_charge':'neutral','amino_hydropathy':4.5,'amino_weight':131.175}
    valine = {'config':('gtg','gtt','gtc','gta'),'name':'valine','letter_code':'V','tri_char':'Val','amino_class':'aliphatic','amino_polarity':'non-polar','amino_charge':'neutral','amino_hydropathy':4.2,'amino_weight':117.148}
    tyrosine = {'config':('tat','tac'),'name':'tyrosine','letter_code':'Y','tri_char':'Tyr','amino_class':'aromatic','amino_polarity':'polar','amino_charge':'neutral','amino_hydropathy':-1.3,'amino_weight':181.191}
    threonine = {'config':('aca','acc','act','acg'),'name':'threonine','letter_code':'T','tri_char':'Thr','amino_class':'hydroxylic','amino_polarity':'polar','amino_charge':'neutral','amino_hydropathy':-0.7,'amino_weight':119.119}
    alanine = {'config':('gcg','gcc','gct','gca'),'name':'alanine','letter_code':'A','tri_char':'Ala','amino_class':'aliphatic','amino_polarity':'non-polar','amino_charge':'neutral','amino_hydropathy':1.8,'amino_weight':89.094}
    asparagine = {'config':('aat','aag'),'name':'asparagine','letter_code':'N','tri_char':'Asn','amino_class':'amide','amino_polarity':'polar','amino_charge':'neutral','amino_hydropathy':-3.5,'amino_weight':132.119}
    lysine = {'config':('aaa','aag'),'name':'lysine','letter_code':'K','tri_char':'Lys','amino_class':'basic','amino_polarity':'basic-polar','amino_charge':'positive','amino_hydropathy':-3.9,'amino_weight':146.189}
    glutamic = {'config':('gaa','gag'),'name':'glutamic','letter_code':'E','tri_char':'Glu','amino_class':'acidic','amino_polarity':'acidic-polar','amino_charge':'negative','amino_hydropathy':-3.5,'amino_weight':147.131}
    aspartic = {'config':('gat','gac'),'name':'aspartic','letter_code':'D','tri_char':'Asp','amino_class':'acidic','amino_polarity':'acidic-polar','amino_charge':'negative','amino_hydropathy':-3.5,'amino_weight':133.104}
    glycine = {'config':('ggg','gga','ggc','ggt'),'name':'glycine','letter_code':'G','tri_char':'Gly','amino_class':'aliphatic','amino_polarity':'non-polar','amino_charge':'neutral','amino_hydropathy':-0.4,'amino_weight':133.104}
    ochre = {'config':('taa'),'name':'Ochre','letter_code':'X','tri_char':'_Ochre','amino_class':'Stop','amino_polarity':None,'amino_charge':None,'amino_hydropathy':None,'amino_weight':0}
    amber = {'config':('tag'),'name':'Amber','letter_code':'B','tri_char':'_Amber','amino_class':'Stop','amino_polarity':None,'amino_charge':None,'amino_hydropathy':None,'amino_weight':0}
    opal = {'config':('tga'),'name':'Opal','letter_code':'Z','tri_char':'_Opal','amino_class':'Stop','amino_polarity':None,'amino_charge':None,'amino_hydropathy':None,'amino_weight':0}
    amino_acids = (methionine,leucine,phenylalanine,serine,cysteine,tryptophan,proline,histidine,glutamine,arginine,isoleucine,valine,tyrosine,threonine,alanine,asparagine,lysine,glutamic,aspartic,glycine,ochre,amber,opal)
    blosum_62 = {'A':{'A':4, 'R':-1,'N':-2,'D':-2,'C':0, 'Q':-1,'E':-1,'G':0, 'H':-2,'I':-1,'L':-1,'K':-1,'M':-1,'F':-2,'P':-1,'S':1, 'T':0, 'W':-3,'Y':-2,'V':0, 'X':0,'B':0,'Z':0},
                 'R':{'A':-1,'R':5, 'N':0, 'D':-2,'C':-3,'Q':1, 'E':0, 'G':-2,'H':0, 'I':-3,'L':-2,'K':2, 'M':-1,'F':-3,'P':-2,'S':-1,'T':-1,'W':-3,'Y':-2,'V':-3,'X':0,'B':0,'Z':0},
                 'N':{'A':-2,'R':0, 'N':6, 'D':1, 'C':-3,'Q':0, 'E':0, 'G':0, 'H':1, 'I':-3,'L':-3,'K':0, 'M':-2,'F':-3,'P':-2,'S':1, 'T':0, 'W':-4,'Y':-2,'V':-3,'X':0,'B':0,'Z':0},
                 'D':{'A':-2,'R':-2,'N':1, 'D':6, 'C':-3,'Q':0, 'E':2, 'G':-1,'H':-1,'I':-3,'L':-4,'K':-1,'M':-3,'F':-3,'P':-1,'S':-1,'T':-1,'W':-2,'Y':-2,'V':-1,'X':0,'B':0,'Z':0},
                 'C':{'A':0, 'R':-3,'N':-3,'D':-3,'C':9, 'Q':-3,'E':-4,'G':-3,'H':-3,'I':-1,'L':-1,'K':-3,'M':-1,'F':-2,'P':-3,'S':-1,'T':-1,'W':-2,'Y':-2,'V':-1,'X':0,'B':0,'Z':0},
                 'Q':{'A':-1,'R':1, 'N':0, 'D':0, 'C':-3,'Q':5, 'E':2, 'G':-2,'H':0, 'I':-3,'L':-2,'K':1, 'M':0, 'F':-3,'P':-1,'S':0, 'T':-1,'W':-2,'Y':-1,'V':-2,'X':0,'B':0,'Z':0},
                 'E':{'A':-1,'R':0, 'N':0, 'D':2, 'C':-4,'Q':2, 'E':5, 'G':-2,'H':0, 'I':-3,'L':-3,'K':1, 'M':-2,'F':-3,'P':-1,'S':0, 'T':-1,'W':-3,'Y':-2,'V':-2,'X':0,'B':0,'Z':0},
                 'G':{'A':0, 'R':-2,'N':0, 'D':-1,'C':-3,'Q':-2,'E':-2,'G':6, 'H':-2,'I':-4,'L':-4,'K':-2,'M':-3,'F':-3,'P':-2,'S':0, 'T':-2,'W':-2,'Y':-3,'V':-3,'X':0,'B':0,'Z':0},
                 'H':{'A':-2,'R':0, 'N':1, 'D':-1,'C':-3,'Q':0, 'E':0, 'G':-2,'H':8, 'I':-3,'L':-3,'K':-1,'M':-2,'F':-1,'P':-2,'S':-1,'T':-2,'W':-2,'Y':2, 'V':-3,'X':0,'B':0,'Z':0},
                 'I':{'A':-1,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-3,'E':-3,'G':-4,'H':-3,'I':4, 'L':2, 'K':-3,'M':1, 'F':0, 'P':-3,'S':-2,'T':-1,'W':-3,'Y':-1,'V':3, 'X':0,'B':0,'Z':0},
                 'L':{'A':-1,'R':-2,'N':-3,'D':-4,'C':-1,'Q':-2,'E':-3,'G':-4,'H':-3,'I':2, 'L':4, 'K':-2,'M':2, 'F':0, 'P':-3,'S':-2,'T':-1,'W':-2,'Y':-1,'V':1, 'X':0,'B':0,'Z':0},
                 'K':{'A':-1,'R':2, 'N':0, 'D':-1,'C':-3,'Q':1, 'E':1, 'G':-2,'H':-1,'I':-3,'L':-2,'K':5, 'M':-1,'F':-3,'P':-1,'S':0, 'T':-1,'W':-3,'Y':-2,'V':-2,'X':0,'B':0,'Z':0},
                 'M':{'A':-1,'R':-1,'N':-2,'D':-3,'C':-1,'Q':0, 'E':-2,'G':-3,'H':-2,'I':1, 'L':2, 'K':-1,'M':5, 'F':0, 'P':-2,'S':-1,'T':-1,'W':-1,'Y':-1,'V':1, 'X':0,'B':0,'Z':0},
                 'F':{'A':-2,'R':-3,'N':-3,'D':-3,'C':-2,'Q':-3,'E':-3,'G':-3,'H':-1,'I':0, 'L':0, 'K':-3,'M':0, 'F':6, 'P':-4,'S':-2,'T':-2,'W':1, 'Y':3, 'V':-1,'X':0,'B':0,'Z':0},
                 'P':{'A':-1,'R':-2,'N':-2,'D':-1,'C':-3,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-3,'L':-3,'K':-1,'M':-2,'F':-4,'P':7, 'S':-1,'T':-1,'W':-4,'Y':-2,'V':-2,'X':0,'B':0,'Z':0},
                 'S':{'A':1, 'R':-1,'N':1, 'D':0, 'C':-1,'Q':0, 'E':0, 'G':0, 'H':-1,'I':-2,'L':-2,'K':0, 'M':-1,'F':-2,'P':-1,'S':4, 'T':1, 'W':-3,'Y':-2,'V':-2,'X':0,'B':0,'Z':0},
                 'T':{'A':0, 'R':-1,'N':0, 'D':-1,'C':-1,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-1,'L':-1,'K':-1,'M':-1,'F':-2,'P':-1,'S':1, 'T':5, 'W':-2,'Y':-2,'V':0, 'X':0,'B':0,'Z':0},
                 'W':{'A':-3,'R':-3,'N':-4,'D':-4,'C':-2,'Q':-2,'E':-3,'G':-2,'H':-2,'I':-3,'L':-2,'K':-3,'M':-1,'F':1, 'P':-4,'S':-3,'T':-2,'W':11,'Y':-2,'V':-3,'X':0,'B':0,'Z':0},
                 'Y':{'A':-2,'R':-2,'N':-2,'D':-3,'C':-2,'Q':-1,'E':-2,'G':-3,'H':2, 'I':-1,'L':-1,'K':-2,'M':-1,'F':3, 'P':-3,'S':-2,'T':-2,'W':-2,'Y':7, 'V':-1,'X':0,'B':0,'Z':0},
                 'V':{'A':0, 'R':-3,'N':-3,'D':-3,'C':-1,'Q':-2,'E':-2,'G':-3,'H':-3,'I':3, 'L':1, 'K':-2,'M':1, 'F':-1,'P':-2,'S':-2,'T':0, 'W':-3,'Y':-1,'V':4, 'X':0,'B':0,'Z':0},       
                 'X':{'A':0,'R':0,'N':0,'D':0,'C':0,'Q':0,'E':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0,'X':0,'B':0,'Z':0},
                 'B':{'A':0,'R':0,'N':0,'D':0,'C':0,'Q':0,'E':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0,'X':0,'B':0,'Z':0},
                 'Z':{'A':0,'R':0,'N':0,'D':0,'C':0,'Q':0,'E':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0,'X':0,'B':0,'Z':0}}    
    def __init__(self):
        for i in range(len(self.amino_acids)-3):
            self.amino_acid.append(self.amino_acids[i]['letter_code'])
class Prediction:
    def __init__(self, possible_range):
        self.possible_range = possible_range
        self.possible_config = []
        self.amino_data = Amino_Data()
        for i in range(possible_range):
            self.possible_config.append([])
    def aa_predict(self, input_seq):
        for i in range(len(input_seq)):
            if input_seq[i] == 'M':
                for i in range(self.possible_range):
                    self.possible_config[i].append('atg')
            elif input_seq[i] == 'L':
                for i in range(6):
                    self.possible_config[i].append(self.amino_data.leucine['config'][i])
                for i in range(self.possible_range-6):
                    self.possible_config[i+6].append(''.join(choices(['tta','ttg','ctt','ctc','cta','ctg'],k=1)))
            elif input_seq[i] == 'F':
                self.possible_config[0].append('ttt')
                self.possible_config[1].append('ttc')
                for i in range(self.possible_range-2):
                    self.possible_config[i+2].append(''.join(choices(['ttt','ttc'],k=1)))
            elif input_seq[i] == 'S':
                for i in range(6):
                    self.possible_config[i].append(self.amino_data.serine['config'][i])
                for i in range(self.possible_range-6):
                    self.possible_config[i+6].append(''.join(choices(['tct','tcc','tca','tcg','agc','agt'],k=1)))
            elif input_seq[i] == 'Y':
                self.possible_config[0].append('tat')
                self.possible_config[1].append('tac')
                for i in range(self.possible_range-2):
                    self.possible_config[i+2].append(''.join(choices(['tat','tac'],k=1)))
            elif input_seq[i] == 'C':
                self.possible_config[0].append('tgt')
                self.possible_config[1].append('tgc')
                for i in range(self.possible_range-2):
                    self.possible_config[i+2].append(''.join(choices(['tgt','tgc'],k=1)))
            elif input_seq[i] == 'W':
                for i in range(self.possible_range):
                    self.possible_config[i].append('tgg')
            elif input_seq[i] == 'P':
                for i in range(4):
                    self.possible_config[i].append(self.amino_data.proline['config'][i])                  
                for i in range(self.possible_range-4):
                    self.possible_config[i+4].append(''.join(choices(['ccc','cca','cct','ccg'],k=1)))
            elif input_seq[i] == 'H':
                self.possible_config[0].append('cat')
                self.possible_config[1].append('cac')
                for i in range(self.possible_range-2):
                    self.possible_config[i+2].append(''.join(choices(['cac','cat'],k=1)))
            elif input_seq[i] == 'Q':
                self.possible_config[0].append('caa')
                self.possible_config[1].append('cag')
                for i in range(self.possible_range-2):
                    self.possible_config[i+2].append(''.join(choices(['caa','cag'],k=1)))
            elif input_seq[i] == 'R':
                for i in range(6):
                    self.possible_config[i].append(self.amino_data.arginine['config'][i])
                for i in range(self.possible_range-6):
                    self.possible_config[i+6].append(''.join(choices(['cgg','cgt','cgc','cga','agg','aga'],k=1)))
            elif input_seq[i] == 'I':
                self.possible_config[0].append('att')
                self.possible_config[1].append('atc')
                self.possible_config[2].append('ata')
                for i in range(self.possible_range-3):
                    self.possible_config[i+3].append(''.join(choices(['att','atc','ata'],k=1)))
            elif input_seq[i] == 'V':
                for i in range(4):
                    self.possible_config[i].append(self.amino_data.valine['config'][i])               
                for i in range(self.possible_range-4):
                    self.possible_config[i+4].append(''.join(choices(['gtg','gtt','gtc','gta'],k=1)))
            elif input_seq[i] == 'T':
                for i in range(4):
                    self.possible_config[i].append(self.amino_data.threonine['config'][i])                 
                for i in range(self.possible_range-4):
                    self.possible_config[i+4].append(''.join(choices(['aca','acc','act','acg'],k=1)))
            elif input_seq[i] == 'A':
                for i in range(4):
                    self.possible_config[i].append(self.amino_data.alanine['config'][i])                
                for i in range(self.possible_range-4):
                    self.possible_config[i+4].append(''.join(choices(['gcg','gcc','gct','gca'],k=1)))
            elif input_seq[i] == 'N':
                self.possible_config[0].append('aat')
                self.possible_config[1].append('aac')
                for i in range(self.possible_range-2):
                    self.possible_config[i+2].append(''.join(choices(['aat','aac'],k=1)))
            elif input_seq[i] == 'K':
                self.possible_config[0].append('aaa')
                self.possible_config[1].append('aag')
                for i in range(self.possible_range-2):
                    self.possible_config[i+2].append(''.join(choices(['aaa','aag'],k=1)))
            elif input_seq[i] == 'E':
                self.possible_config[0].append('gaa')
                self.possible_config[1].append('gag')
                for i in range(self.possible_range-2):
                    self.possible_config[i+2].append(''.join(choices(['gaa','gag'],k=1)))
            elif input_seq[i] == 'D':
                self.possible_config[0].append('gat')
                self.possible_config[1].append('gac')
                for i in range(self.possible_range-2):
                    self.possible_config[i+2].append(''.join(choices(['gat','gac'],k=1)))
            elif input_seq[i] == 'G':
                for i in range(4):
                    self.possible_config[i].append(self.amino_data.glycine['config'][i])
                for i in range(self.possible_range-4):
                    self.possible_config[i+4].append(''.join(choices(['ggg','gga','ggc','ggt'],k=1)))
            elif input_seq[i] == 'X':
                for i in range(self.possible_range):
                    self.possible_config[i].append('taa')
            elif input_seq[i] == 'B':
                for i in range(self.possible_range):
                    self.possible_config[i].append('tag')
            elif input_seq[i] == 'Z':
                for i in range(self.possible_range):
                    self.possible_config[i].append('tga')
class Translation:    
    dna_data = Nucleotide_Data()
    amino_data = Amino_Data()
    def translate(self, input_seq):
        self.dna_typed = []
        self.rna_typed = []
        self.dna_match = []
        for i in range(len(input_seq)):
            for f in range(len(self.dna_data.nucleotide)):
                if input_seq[i] in self.dna_data.nucleotide[f]['char']:
                    self.dna_typed.append(self.dna_data.nucleotide[f]['dna'])
                    self.rna_typed.append(self.dna_data.nucleotide[f]['rna'])
                    self.dna_match.append(self.dna_data.nucleotide[f]['match'])
        self.direct_start_index = []
        self.reversed_start_index = []
        self.direct_sequences = []
        self.possible_chain = {}
        self.stopped_chains = {}
        self.stopped_chains_list = []
        self.translated_chains_ammount = []
        self.translated_chains_letter = []
        self.translated_chains_chars = []
        self.translated_chains_weight = []
        self.translated_chains_total_weight = []        
        self.translated_chains_blosum_score = []
        self.translated_chains_blosum_total_score = []
        self.translated_chains_hydropathy = []
        self.translated_chains_class = []
        self.translated_chains_polarity = []
        self.translated_chains_charge = []
        self.translated_chains_name = []
        self.chain_sequence = ''.join(self.dna_typed)
        self.dna_typed.reverse()
        self.reversed_sequence = ''.join(self.dna_typed)
        self.dna_typed.reverse()
        for i in range(len(self.chain_sequence)):
            if self.chain_sequence[i:i+3] == self.amino_data.methionine['config']:
                self.direct_start_index.append(i)
        for i in range(len(self.reversed_sequence)):
            if self.reversed_sequence[i:i+3] == self.amino_data.methionine['config']:
                self.reversed_start_index.append(i)
        for i in range(len(self.direct_start_index)):
            index_var = self.direct_start_index[i]
            self.direct_sequences.append(self.chain_sequence[index_var:])
        for i in range(len(self.reversed_start_index)):
            index_var = self.reversed_start_index[i]
            self.direct_sequences.append(self.reversed_sequence[index_var:])
        for i in range(len(self.direct_sequences)):
            chain_names = "chain{}".format(i)
            for f in range(len(self.direct_sequences[i])):
                self.possible_chain[chain_names] = [self.direct_sequences[i][f:f+3] for f in range(0, len(self.direct_sequences[i]), 3)]
            for f in range(len(self.possible_chain[chain_names])):
                codon = self.possible_chain[chain_names][f]
                if codon == self.amino_data.opal['config'] or codon == self.amino_data.amber['config'] or codon == self.amino_data.ochre['config']:
                    stop_index = self.possible_chain[chain_names].index(codon)+1
                    self.stopped_chains[chain_names] = self.possible_chain[chain_names][0:stop_index]
        for i in self.stopped_chains:
            self.stopped_chains_list.append(self.stopped_chains[i])
        for i in range(len(self.stopped_chains_list)):
            self.translated_chains_ammount.append(len(self.stopped_chains_list[i]))
            self.translated_chains_letter.append([])
            self.translated_chains_chars.append([])
            self.translated_chains_weight.append([])
            self.translated_chains_total_weight += [0]
            self.translated_chains_blosum_score.append([])
            self.translated_chains_blosum_total_score += [0]
            self.translated_chains_hydropathy.append([])
            self.translated_chains_class.append([])
            self.translated_chains_polarity.append([])
            self.translated_chains_charge.append([])
            self.translated_chains_name.append([])
            for f in range(len(self.stopped_chains_list[i])):
                for k in range(len(self.amino_data.amino_acids)):
                    if self.stopped_chains_list[i][f] in self.amino_data.amino_acids[k]['config']:
                        self.translated_chains_letter[i].append(self.amino_data.amino_acids[k]['letter_code'])
                        self.translated_chains_chars[i].append(self.amino_data.amino_acids[k]['tri_char'])
                        self.translated_chains_weight[i].append(self.amino_data.amino_acids[k]['amino_weight'])
                        self.translated_chains_hydropathy[i].append(self.amino_data.amino_acids[k]['amino_hydropathy'])
                        self.translated_chains_class[i].append(self.amino_data.amino_acids[k]['amino_class'])
                        self.translated_chains_polarity[i].append(self.amino_data.amino_acids[k]['amino_polarity'])
                        self.translated_chains_charge[i].append(self.amino_data.amino_acids[k]['amino_charge'])
                        self.translated_chains_name[i].append(self.amino_data.amino_acids[k]['name'])
            for f in self.translated_chains_weight[i]:
                self.translated_chains_total_weight[i] = sum(self.translated_chains_weight[i])
            for f in range(len(self.translated_chains_letter[i])-1):
                self.translated_chains_blosum_score[i].append([])
                first_letter_var = self.translated_chains_letter[i][f]
                second_letter_var = self.translated_chains_letter[i][f+1]
                blosum_score_var = self.amino_data.blosum_62[first_letter_var][second_letter_var]
                self.translated_chains_blosum_score[i][f] = blosum_score_var    
            for f in self.translated_chains_blosum_score[i]:
                self.translated_chains_blosum_total_score[i] = sum(self.translated_chains_blosum_score[i])
class Tk_Main(Tk):
    def __init__(self):
        super().__init__()
        self.title("Biochem tool: DNA Prediction, transcription, translation, and Amino acid evaluation.")
        self.canvas = Canvas(self)
        self.frame = Frame(self.canvas)
        self.draw_index_frame()
    def draw_index_frame(self):
        self.input_entry = StringVar()
        self.entry_var = Entry(self, width=50, text=self.input_entry)
        self.r_dot = StringVar()
        self.entry_var.grid(column=0, row=0)
        self.predict = Radiobutton(self, text='DNA configuration prediction from peptide sequence', command=self.select_predict,variable=self.r_dot, value='predict')
        self.predict.grid(column=0, row=1)
        self.translate = Radiobutton(self, text='DNA/mRNA, translation, and peptide stat', command=self.select_eval,variable=self.r_dot, value='translate')
        self.translate.grid(column=0, row=2)
        self.entry_var.focus()
        self.eval('tk::PlaceWindow %s center' % self.winfo_pathname(self.winfo_id()))
        self.eval_selected = 0
        self.predict_selected = 0
        self.select_eval()
        self.r_dot.set('translate') 
    def remove_first_frame(self):
        Entry.destroy(self.entry_var)
        Radiobutton.destroy(self.predict)
        Radiobutton.destroy(self.translate)
        if self.eval_selected == 1:
            self.remove_eval_frame()
            self.eval_selected = 0
        if self.predict_selected == 1:
            self.remove_predict_frame()
            self.predict_selected = 0    
    def run_analysis(self, *args):    
        if self.r_dot.get() == 'predict':
            input_seq = self.input_entry.get()
            predict_range = int(self.p_dot.get()) + 6
            self.amino_predict = Prediction(predict_range)
            self.amino_predict.aa_predict(input_seq)
            self.remove_first_frame()
            self.draw_prediction_frame()
            self.predict_selected = 0
        elif self.r_dot.get() == 'translate':   
            input_seq = self.input_entry.get()
            self.amino_translate = Translation()
            self.amino_translate.translate(input_seq)
            self.remove_first_frame()
            self.draw_evaluation_frame()
            self.eval_selected = 0
        self.eval('tk::PlaceWindow %s center' % self.winfo_pathname(self.winfo_id()))
    def select_eval(self):
        self.eval_selected += 1
        if self.predict_selected >= 1: 
            self.remove_predict_frame()
            self.predict_selected = 0
        if self.eval_selected == 1:
            self.check_dna_var = StringVar()
            self.check_dna = Checkbutton(self, text='show DNA', variable=self.check_dna_var, onvalue='true', offvalue='false')
            self.check_dna.grid()
            self.check_match_var = StringVar()
            self.check_match = Checkbutton(self, text='show DNA Matching codon', variable=self.check_match_var, onvalue='true', offvalue='false')
            self.check_match.grid()
            self.check_rna_var = StringVar()
            self.check_rna = Checkbutton(self, text='show mRNA', variable=self.check_rna_var, onvalue='true', offvalue='false')
            self.check_rna.grid()
            self.check_chars_var = StringVar()
            self.check_chars = Checkbutton(self, text='show 3 letter code', variable=self.check_chars_var, onvalue='true', offvalue='false')
            self.check_chars.grid()
            self.check_weight_distr_var = StringVar()
            self.check_weight_distr = Checkbutton(self, text='show weight distribution', variable=self.check_weight_distr_var, onvalue='true', offvalue='false')
            self.check_weight_distr.grid()
            self.check_blosum_distr_var = StringVar()
            self.check_blosum_distr = Checkbutton(self, text='show individual BLOSUM62 score', variable=self.check_blosum_distr_var, onvalue='true', offvalue='false')
            self.check_blosum_distr.grid()
            self.check_polarity_distr_var = StringVar()
            self.check_polarity_distr = Checkbutton(self, text='show polarity distribution', variable=self.check_polarity_distr_var, onvalue='true', offvalue='false')
            self.check_polarity_distr.grid()
            self.check_hydropathy_distr_var = StringVar()
            self.check_hydropathy_distr = Checkbutton(self, text='show hydropathy distribution', variable=self.check_hydropathy_distr_var, onvalue='true', offvalue='false')
            self.check_hydropathy_distr.grid()
            self.check_class_distr_var = StringVar()
            self.check_class_distr = Checkbutton(self, text='show class distribution', variable=self.check_class_distr_var, onvalue='true', offvalue='false')
            self.check_class_distr.grid()
            self.check_charge_distr_var = StringVar()
            self.check_charge_distr = Checkbutton(self, text='show charge distribution', variable=self.check_charge_distr_var, onvalue='true', offvalue='false')
            self.check_charge_distr.grid()
            self.eval_btn = Button(self, text='Run translation', command=self.run_analysis)
            self.eval_btn.grid()
            self.bind('<Return>', self.run_analysis)
    def draw_evaluation_frame(self):
        self.transform_chain = str.maketrans(",'", "  ")
        self.transform_data = str.maketrans("'", " ")
        self.chain_ammount_config_lbl = Label(self, text="Chain ammount:")
        self.chain_ammount_config_lbl.grid()
        self.chain_ammount_config_var_txt = len(self.amino_translate.stopped_chains_list)
        self.chain_ammount_config_var_lbl = Label(self, text=self.chain_ammount_config_var_txt)
        self.chain_ammount_config_var_lbl.grid()
        self.amino_ammount_config_lbl = Label(self, text="Amino ammount:")
        self.amino_ammount_config_lbl.grid()
        self.amino_ammount_config_var_txt = self.amino_translate.translated_chains_ammount
        self.amino_ammount_config_var_lbl = Label(self, text=self.amino_ammount_config_var_txt)
        self.amino_ammount_config_var_lbl.grid()
        self.total_mass_config_lbl = Label(self, text="Total mass:")
        self.total_mass_config_lbl.grid()
        self.total_mass_config_var_txt = self.amino_translate.translated_chains_total_weight
        self.total_mass_config_var_lbl = Label(self, text=self.total_mass_config_var_txt)
        self.total_mass_config_var_lbl.grid()
        self.total_blosum_config_lbl = Label(self, text="Total BLOSUM62 score:")
        self.total_blosum_config_lbl.grid()
        self.total_blosum_config_var_txt = self.amino_translate.translated_chains_blosum_total_score
        self.total_blosum_config_var_lbl = Label(self, text=self.total_blosum_config_var_txt)
        self.total_blosum_config_var_lbl.grid()
        self.amino_letter_config_lbl = Label(self, text="Amino letter:")
        self.amino_letter_config_lbl.grid()
        self.amino_letter_config_var_txt_r = str(self.amino_translate.translated_chains_letter)
        self.amino_letter_config_var_tr = str.translate(self.amino_letter_config_var_txt_r, self.transform_chain)
        self.amino_letter_config_var_txt = self.amino_letter_config_var_tr.replace(" ", "")
        self.amino_letter_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
        self.amino_letter_config_var_lbl.insert(float(), self.amino_letter_config_var_txt)
        self.amino_letter_config_var_lbl.grid()
        if self.check_dna_var.get() == 'true':
            self.dna_typed_lbl = Label(self, text="DNA:")
            self.dna_typed_lbl.grid()
            self.dna_typed_var_txt_r = str(self.amino_translate.dna_typed)
            self.dna_typed_var_tr = str.translate(self.dna_typed_var_txt_r, self.transform_chain)
            self.dna_typed_var_txt = self.dna_typed_var_tr.replace(" ", "")
            self.dna_typed_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
            self.dna_typed_var_lbl.insert(float(), self.dna_typed_var_txt)
            self.dna_typed_var_lbl.grid()
        if self.check_match_var.get() == 'true':
            self.match_config_lbl = Label(self, text="DNA Match:")
            self.match_config_lbl.grid()
            self.match_config_var_txt_r = str(self.amino_translate.dna_match)
            self.match_config_var_tr = str.translate(self.match_config_var_txt_r, self.transform_chain)
            self.match_config_var_txt = self.match_config_var_tr.replace(" ", "")
            self.match_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
            self.match_config_var_lbl.insert(float(), self.match_config_var_txt)
            self.match_config_var_lbl.grid()
        if self.check_rna_var.get() == 'true':
            self.rna_config_lbl = Label(self, text="mRNA:")
            self.rna_config_lbl.grid()
            self.rna_config_var_txt_r = str(self.amino_translate.rna_typed)
            self.rna_config_var_tr = str.translate(self.rna_config_var_txt_r, self.transform_chain)
            self.rna_config_var_txt = self.rna_config_var_tr.replace(" ", "")
            self.rna_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
            self.rna_config_var_lbl.insert(float(), self.rna_config_var_txt)
            self.rna_config_var_lbl.grid()
        if self.check_chars_var.get() == 'true':
            self.amino_chars_config_lbl = Label(self, text="Amino chars:")
            self.amino_chars_config_lbl.grid()
            self.amino_chars_config_var_txt_r = str(self.amino_translate.translated_chains_chars)
            self.amino_chars_config_var_tr = str.translate(self.amino_chars_config_var_txt_r, self.transform_chain)
            self.amino_chars_config_var_txt = self.amino_chars_config_var_tr.replace(" ", "")
            self.amino_chars_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
            self.amino_chars_config_var_lbl.insert(float(), self.amino_chars_config_var_txt)
            self.amino_chars_config_var_lbl.grid()
        if self.check_weight_distr_var.get() == 'true':
            self.mass_distr_config_lbl = Label(self, text="Weight distribution:")
            self.mass_distr_config_lbl.grid()
            self.mass_distr_config_var_txt_r = str(self.amino_translate.translated_chains_weight)
            self.mass_distr_config_var_tr = str.translate(self.mass_distr_config_var_txt_r, self.transform_data)
            self.mass_distr_config_var_txt = self.mass_distr_config_var_tr.replace(" ", "")
            self.mass_distr_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
            self.mass_distr_config_var_lbl.insert(float(), self.mass_distr_config_var_txt)
            self.mass_distr_config_var_lbl.grid()
        if self.check_blosum_distr_var.get() == 'true':
            self.blosum_score_config_lbl = Label(self, text="Individual BLOSUM62 score:")
            self.blosum_score_config_lbl.grid()
            self.blosum_score_config_var_txt_r = str(self.amino_translate.translated_chains_blosum_score)
            self.blosum_score_config_var_tr = str.translate(self.blosum_score_config_var_txt_r, self.transform_data)
            self.blosum_score_config_var_txt = self.blosum_score_config_var_tr.replace(" ", "")
            self.blosum_score_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
            self.blosum_score_config_var_lbl.insert(float(), self.blosum_score_config_var_txt)
            self.blosum_score_config_var_lbl.grid()
        if self.check_polarity_distr_var.get() == 'true':
            self.distr_polarity_config_lbl = Label(self, text="Distributed polarity:")
            self.distr_polarity_config_lbl.grid()
            self.distr_polarity_config_var_txt_r = str(self.amino_translate.translated_chains_polarity)
            self.distr_polarity_config_var_tr = str.translate(self.distr_polarity_config_var_txt_r, self.transform_data)
            self.distr_polarity_config_var_txt = self.distr_polarity_config_var_tr.replace(" ", "")
            self.distr_polarity_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
            self.distr_polarity_config_var_lbl.insert(float(), self.distr_polarity_config_var_txt)
            self.distr_polarity_config_var_lbl.grid()
        if self.check_hydropathy_distr_var.get() == 'true':
            self.distr_hydropathy_config_lbl = Label(self, text="Distributed hydropathy:")
            self.distr_hydropathy_config_lbl.grid()
            self.distr_hydropathy_config_var_txt_r = str(self.amino_translate.translated_chains_hydropathy)
            self.distr_hydropathy_config_var_tr = str.translate(self.distr_hydropathy_config_var_txt_r, self.transform_data)
            self.distr_hydropathy_config_var_txt = self.distr_hydropathy_config_var_tr.replace(" ", "")
            self.distr_hydropathy_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
            self.distr_hydropathy_config_var_lbl.insert(float(), self.distr_hydropathy_config_var_txt)
            self.distr_hydropathy_config_var_lbl.grid()
        if self.check_class_distr_var.get() == 'true':
            self.distr_class_config_lbl = Label(self, text="Distributed class:")
            self.distr_class_config_lbl.grid()
            self.distr_class_config_var_txt_r = str(self.amino_translate.translated_chains_class)
            self.distr_class_config_var_tr = str.translate(self.distr_class_config_var_txt_r, self.transform_data)
            self.distr_class_config_var_txt = self.distr_class_config_var_tr.replace(" ", "")
            self.distr_class_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
            self.distr_class_config_var_lbl.insert(float(), self.distr_class_config_var_txt)
            self.distr_class_config_var_lbl.grid()
        if self.check_charge_distr_var.get() == 'true':
            self.distr_charge_config_lbl = Label(self, text="Distributed charge:")
            self.distr_charge_config_lbl.grid()
            self.distr_charge_config_var_txt_r = str(self.amino_translate.translated_chains_charge)
            self.distr_charge_config_var_tr = str.translate(self.distr_charge_config_var_txt_r, self.transform_data)
            self.distr_charge_config_var_txt = self.distr_charge_config_var_tr.replace(" ", "")
            self.distr_charge_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
            self.distr_charge_config_var_lbl.insert(float(), self.distr_charge_config_var_txt)
            self.distr_charge_config_var_lbl.grid()
        self.return_btn = Button(self, text='Go Back', command=self.remove_translation_frame)
        self.return_btn.grid()
        self.bind('<Return>', self.remove_translation_frame)
    def remove_eval_frame(self):
        Checkbutton.destroy(self.check_dna)
        Checkbutton.destroy(self.check_match)
        Checkbutton.destroy(self.check_rna)
        Checkbutton.destroy(self.check_chars)
        Checkbutton.destroy(self.check_weight_distr)
        Checkbutton.destroy(self.check_blosum_distr)
        Checkbutton.destroy(self.check_polarity_distr)
        Checkbutton.destroy(self.check_hydropathy_distr)
        Checkbutton.destroy(self.check_class_distr)
        Checkbutton.destroy(self.check_charge_distr)
        Button.destroy(self.eval_btn)
    def remove_translation_frame(self, *args):
        Label.destroy(self.chain_ammount_config_lbl)
        Text.destroy(self.chain_ammount_config_var_lbl)
        Label.destroy(self.amino_ammount_config_lbl)
        Text.destroy(self.amino_ammount_config_var_lbl)
        Label.destroy(self.total_mass_config_lbl)
        Text.destroy(self.total_mass_config_var_lbl)
        Label.destroy(self.total_blosum_config_lbl)
        Text.destroy(self.total_blosum_config_var_lbl)
        Label.destroy(self.amino_letter_config_lbl)
        Text.destroy(self.amino_letter_config_var_lbl)
        if self.check_dna_var.get() == 'true':
            Label.destroy(self.dna_typed_lbl)
            Text.destroy(self.dna_typed_var_lbl)
        if self.check_match_var.get() == 'true':
            Label.destroy(self.match_config_lbl)
            Text.destroy(self.match_config_var_lbl)
        if self.check_rna_var.get() == 'true':
            Label.destroy(self.rna_config_lbl)
            Text.destroy(self.rna_config_var_lbl)
        if self.check_chars_var.get() == 'true':
            Label.destroy(self.amino_chars_config_lbl)
            Text.destroy(self.amino_chars_config_var_lbl)
        if self.check_weight_distr_var.get() == 'true':
            Label.destroy(self.mass_distr_config_lbl)
            Text.destroy(self.mass_distr_config_var_lbl)
        if self.check_blosum_distr_var.get() == 'true':
            Label.destroy(self.blosum_score_config_lbl)
            Text.destroy(self.blosum_score_config_var_lbl)
        if self.check_polarity_distr_var.get() == 'true':
            Label.destroy(self.distr_polarity_config_lbl)
            Text.destroy(self.distr_polarity_config_var_lbl)
        if self.check_hydropathy_distr_var.get() == 'true':
            Label.destroy(self.distr_hydropathy_config_lbl)
            Text.destroy(self.distr_hydropathy_config_var_lbl)
        if self.check_class_distr_var.get() == 'true':
            Label.destroy(self.distr_class_config_lbl)
            Text.destroy(self.distr_class_config_var_lbl)
        if self.check_charge_distr_var.get() == 'true':
            Label.destroy(self.distr_charge_config_lbl)
            Text.destroy(self.distr_charge_config_var_lbl)
        Button.destroy(self.return_btn)
        self.draw_index_frame()
    def select_predict(self):
        self.predict_selected += 1
        if self.eval_selected >= 1:
            self.remove_eval_frame()
            self.eval_selected = 0
        if self.predict_selected == 1:
            self.p_dot = StringVar()
            self.six_predict = Radiobutton(self, text='6 possible linear settings', variable=self.p_dot, value='0')
            self.six_predict.grid()
            self.seven_predict = Radiobutton(self, text='6 linear settings + 1 random', variable=self.p_dot, value='1')
            self.seven_predict.grid()
            self.eight_predict = Radiobutton(self, text='6 linear settings + 2 random', variable=self.p_dot, value='2')
            self.eight_predict.grid()
            self.nine_predict = Radiobutton(self, text='6 linear settings + 3 random', variable=self.p_dot, value='3')
            self.nine_predict.grid()
            self.ten_predict = Radiobutton(self, text='6 linear settings + 4 random', variable=self.p_dot, value='4')
            self.ten_predict.grid()
            self.predict_btn = Button(self, text='Run prediction', command=self.run_analysis)
            self.p_dot.set('0')
            self.predict_btn.grid()   
            self.bind('<Return>', self.run_analysis)
    def remove_predict_frame(self):
        Radiobutton.destroy(self.six_predict)
        Radiobutton.destroy(self.seven_predict)
        Radiobutton.destroy(self.eight_predict)
        Radiobutton.destroy(self.nine_predict)
        Radiobutton.destroy(self.ten_predict)
        Button.destroy(self.predict_btn)
    def draw_prediction_frame(self):
        self.transform_chain = str.maketrans(",'", "  ")       
        self.first_config_lbl = Label(self, text="Linear config 1:")
        self.first_config_lbl.grid()
        self.first_config_txt_r = str(self.amino_predict.possible_config[0])
        self.first_config_tr = str.translate(self.first_config_txt_r, self.transform_chain)
        self.first_config_txt = self.first_config_tr.replace(" ", "")
        self.first_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
        self.first_config_var_lbl.insert(float(), self.first_config_txt)
        self.first_config_var_lbl.grid()
        self.second_config_lbl = Label(self, text="Linear config 2:")
        self.second_config_lbl.grid()
        self.second_config_txt_r = str(self.amino_predict.possible_config[1])
        self.second_config_tr = str.translate(self.second_config_txt_r, self.transform_chain)
        self.second_config_txt = self.second_config_tr.replace(" ", "")
        self.second_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
        self.second_config_var_lbl.insert(float(), self.second_config_txt)
        self.second_config_var_lbl.grid()
        self.third_config_lbl = Label(self, text="Linear config 3:")
        self.third_config_lbl.grid()
        self.third_config_txt_r = str(self.amino_predict.possible_config[2])
        self.third_config_tr = str.translate(self.third_config_txt_r, self.transform_chain)
        self.third_config_txt = self.third_config_tr.replace(" ", "")
        self.third_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
        self.third_config_var_lbl.insert(float(), self.third_config_txt)
        self.third_config_var_lbl.grid()
        self.fourth_config_lbl = Label(self, text="Linear config 4:")
        self.fourth_config_lbl.grid()
        self.fourth_config_txt_r = str(self.amino_predict.possible_config[3])
        self.fourth_config_tr = str.translate(self.fourth_config_txt_r, self.transform_chain)
        self.fourth_config_txt = self.fourth_config_tr.replace(" ", "")
        self.fourth_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
        self.fourth_config_var_lbl.insert(float(), self.fourth_config_txt)
        self.fourth_config_var_lbl.grid()
        self.fifth_config_lbl = Label(self, text="Linear config 5:")
        self.fifth_config_lbl.grid()
        self.fifth_config_txt_r = str(self.amino_predict.possible_config[4])
        self.fifth_config_tr = str.translate(self.fifth_config_txt_r, self.transform_chain)
        self.fifth_config_txt = self.fifth_config_tr.replace(" ", "")
        self.fifth_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
        self.fifth_config_var_lbl.insert(float(), self.fifth_config_txt)
        self.fifth_config_var_lbl.grid()
        self.sixth_config_lbl = Label(self, text="Linear config 6:")
        self.sixth_config_lbl.grid()
        self.sixth_config_txt_r = str(self.amino_predict.possible_config[5])
        self.sixth_config_tr = str.translate(self.sixth_config_txt_r, self.transform_chain)
        self.sixth_config_txt = self.sixth_config_tr.replace(" ", "")
        self.sixth_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
        self.sixth_config_var_lbl.insert(float(), self.sixth_config_txt)
        self.sixth_config_var_lbl.grid()
        if int(self.p_dot.get()) == 1:
            self.declare_rand_conf_1()
        if int(self.p_dot.get()) == 2:
            self.declare_rand_conf_1()
            self.declare_rand_conf_2()
        if int(self.p_dot.get()) == 3:
            self.declare_rand_conf_1()
            self.declare_rand_conf_2()
            self.declare_rand_conf_3()
        if int(self.p_dot.get()) == 4:
            self.declare_rand_conf_1()
            self.declare_rand_conf_2()
            self.declare_rand_conf_3()
            self.declare_rand_conf_4()
        self.return_btn = Button(self, text='Go back', command=self.remove_configuration_frame)
        self.return_btn.grid()
        self.bind('<Return>', self.remove_configuration_frame)
    def declare_rand_conf_1(self):
        self.eight_config_lbl = Label(self, text="Randomised config 7:")
        self.eight_config_lbl.grid()
        self.eight_config_txt_r = str(self.amino_predict.possible_config[6])
        self.eight_config_tr = str.translate(self.eight_config_txt_r, self.transform_chain)
        self.eight_config_txt = self.eight_config_tr.replace(" ", "")
        self.eight_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
        self.eight_config_var_lbl.insert(float(), self.eight_config_txt)
        self.eight_config_var_lbl.grid()
    def declare_rand_conf_2(self):
        self.nine_config_lbl = Label(self, text="Randomised config 8:")
        self.nine_config_lbl.grid()
        self.nine_config_txt_r = str(self.amino_predict.possible_config[7])
        self.nine_config_tr = str.translate(self.nine_config_txt_r, self.transform_chain)
        self.nine_config_txt = self.nine_config_tr.replace(" ", "")
        self.nine_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
        self.nine_config_var_lbl.insert(float(), self.nine_config_txt)
        self.nine_config_var_lbl.grid()
    def declare_rand_conf_3(self):
        self.ten_config_lbl = Label(self, text="Randomised config 9")
        self.ten_config_lbl.grid()
        self.ten_config_txt_r = str(self.amino_predict.possible_config[8])
        self.ten_config_tr = str.translate(self.ten_config_txt_r, self.transform_chain)
        self.ten_config_txt = self.ten_config_tr.replace(" ", "")
        self.ten_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
        self.ten_config_var_lbl.insert(float(), self.ten_config_txt)
        self.ten_config_var_lbl.grid()
    def declare_rand_conf_4(self):
        self.eleven_config_lbl = Label(self, text="Randomised config 10:")
        self.eleven_config_lbl.grid()
        self.eleven_config_txt_r = str(self.amino_predict.possible_config[9])
        self.eleven_config_tr = str.translate(self.eleven_config_txt_r, self.transform_chain)
        self.eleven_config_txt = self.eleven_config_tr.replace(" ", "")
        self.eleven_config_var_lbl = Text(self, height=1, wrap="none", font=("Helvetica", 14))
        self.eleven_config_var_lbl.insert(float(), self.eleven_config_txt)
        self.eleven_config_var_lbl.grid()
    def remove_configuration_frame(self, *args):
        Label.destroy(self.first_config_lbl)
        Text.destroy(self.first_config_var_lbl)
        Label.destroy(self.second_config_lbl)
        Text.destroy(self.second_config_var_lbl)
        Label.destroy(self.third_config_lbl)
        Text.destroy(self.third_config_var_lbl)
        Label.destroy(self.fourth_config_lbl)
        Text.destroy(self.fourth_config_var_lbl)
        Label.destroy(self.fifth_config_lbl)
        Text.destroy(self.fifth_config_var_lbl)
        Label.destroy(self.sixth_config_lbl)
        Text.destroy(self.sixth_config_var_lbl)
        if int(self.p_dot.get()) == 1:
            Label.destroy(self.eight_config_lbl)
            Text.destroy(self.eight_config_var_lbl)
        if int(self.p_dot.get()) == 2:
            Label.destroy(self.eight_config_lbl)
            Text.destroy(self.eight_config_var_lbl)
            Label.destroy(self.nine_config_lbl)
            Text.destroy(self.nine_config_var_lbl)
        if int(self.p_dot.get()) == 3:
            Label.destroy(self.eight_config_lbl)
            Text.destroy(self.eight_config_var_lbl)
            Label.destroy(self.nine_config_lbl)
            Text.destroy(self.nine_config_var_lbl)
            Label.destroy(self.ten_config_lbl)
            Text.destroy(self.ten_config_var_lbl)
        if int(self.p_dot.get()) == 4:
            Label.destroy(self.eight_config_lbl)
            Text.destroy(self.eight_config_var_lbl)
            Label.destroy(self.nine_config_lbl)
            Text.destroy(self.nine_config_var_lbl)
            Label.destroy(self.ten_config_lbl)
            Text.destroy(self.ten_config_var_lbl)
            Label.destroy(self.eleven_config_lbl)
            Text.destroy(self.eleven_config_var_lbl)
        Button.destroy(self.return_btn)
        self.draw_index_frame()
tk_init = Tk_Main()
tk_init.mainloop()
