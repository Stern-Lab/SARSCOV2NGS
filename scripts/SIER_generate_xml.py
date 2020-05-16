#!/usr/bin/env python3
from Bio import SeqIO
import argparse


# TODO comment default values
# TODO SET INIT VALUES
parser = argparse.ArgumentParser()
parser.add_argument('--final_aln', help='fasta alignment',
                    default='')
parser.add_argument('--xmlTemplate', help='template xml file',
                    default='')
parser.add_argument('--T0', help='begining of time series', default=2019.7)
parser.add_argument('--sirT0', help='begining of epi dynamics',
                    default=2020.085)       # February 1
parser.add_argument('--s', help='number of susceptible individuals',
                    default=9000000)       # February 1
parser.add_argument('--latPeriod',
                    help='duration of exposed (latent) period (days)',
                    default=3)       # 10.1126/science.abb6105
parser.add_argument('--infPeriod',
                    help='duration of infecitous period (days)', default=4)
parser.add_argument('--pH',
                    help='proportion of I_high individuals', default=0.2)
parser.add_argument('--P',
                    help='proportion of infection caused by I_high individuals',
                    default=0.85)   # JLS superspreading paper, TODO find doi
parser.add_argument('--R0M',
                    help='desired mean for lognormal R0 prior', default=1.5)
parser.add_argument('--R0S',
                    help='desired s for lognormal R0 prior', default=0.5)
parser.add_argument('--R0Change',
                    help='date on which R0 changes by factor a', default='NA')    # default is no R0 change
parser.add_argument('--exogGrM',
                    help='exog growth rate lognormal mean', default=3.6)
parser.add_argument('--exogGrS',
                    help='exog growth rate lognormal standard deviation',
                    default=1)
parser.add_argument('--importRateM',
                    help='mean of import rate prior (exponential)', default=10)
parser.add_argument('--importRateU',
                    help='upper limit of import rate', default=10)
parser.add_argument('--importChangeDate',
                    help='date when importation changes', default=2050) # default no importation change
parser.add_argument('--importChange',
                    help='factor by which import rate changes', default=1)
parser.add_argument('--initTree',
                    help='starting tree',
                    default='')
parser.add_argument('--out_name',
                    help='output file name', default='formatted')
args = parser.parse_args()
args.base_dir = '/'.join(args.final_aln.split('/')[0:-1])
args.base_path = args.base_dir+'/'+args.out_name


def generate_sequence_block(sequence_list):
    sequence_block = ''
    for item in sequence_list:
        sequence_block += '<sequence id="seq_{0}" spec="Sequence" taxon="{0}" \
        totalcount="4" value="{1}" />\n'.format(item.id, item.seq)
    return(sequence_block)


def generate_time_block(sequence_list):
    time_block = '"'
    for item in sequence_list:
        time_block += '{0}={1},'.format(item.id, item.id.split('_')[-2])
    # removes trailing ','
    time_block = ''.join(list(time_block)[0:-2])
    time_block += '"'
    return(time_block)


def generate_xml(args):
    masked_fasta = list(SeqIO.parse(args.finalAln, 'fasta'))
    sequence_block = generate_sequence_block(masked_fasta)
    time_block = generate_time_block(masked_fasta)
    tau = (1-args.pH)*(1/args.pH)/(1/args.P - 1)
    with open(args.finalTree, 'r') as infile:
        initTree = infile.read()
    with open(args.xmlTemplate, 'r') as infile:
        template = infile.read()
    template = template.replace('SEIR_TEMPLATE', args.out_name)
    template = template.replace('<!-- SEQUENCE_BLOCK -->', sequence_block)
    template = template.replace('<!-- TIME_BLOCK -->', time_block)
    template = template.replace('<!-- TREE_BLOCK -->', '"' + initTree + '"')
    template = template.replace('<!-- T0 -->', str(args.T0))
    template = template.replace('<!-- SIR_T0 -->', str(args.sirT0))
    template = template.replace('<!-- S -->', str(args.s))
    template = template.replace('<!-- LATENT_PERIOD -->',
                                str(365/args.latPeriod))
    template = template.replace('<!-- INF_PERIOD -->', str(365/args.infPeriod))
    template = template.replace('<!-- PH -->', str(args.pH))
    template = template.replace('<!-- TAU -->', str(tau))
    template = template.replace('<!-- R0M -->', str(args.R0M))
    template = template.replace('<!-- R0S -->', str(args.R0S))
    template = template.replace('<!-- EXOG_PERIOD -->',
                                str(365/(args.infPeriod+args.latPeriod)))
    template = template.replace('<!-- EXOG_INIT -->',
                                str(args.exogInit))
    template = template.replace('<!-- EXOG_INIT_M -->',
                                str(args.exogInitM))
    template = template.replace('<!-- EXOG_GROWTH_RATE_M -->',
                                str(args.exogGrM))
    template = template.replace('<!-- EXOG_GROWTH_RATE_S -->',
                                str(args.exogGrS))
    template = template.replace('<!-- IMPORT_RATE_M -->',
                                str(args.importRateM))
    template = template.replace('<!-- IMPORT_RATE_U -->',
                                str(args.importRateU))
    template = template.replace('<!-- IMPORT_CHANGE_DATE -->',
                                str(args.importChangeDate))
    template = template.replace('<!-- IMPORT_CHANGE -->',
                                str(args.importChange))
    if args.R0Change == 'NA':
        # todo make this...better?
        R0Change_operator = '<operator id="seir.a.operator.t:{0}" \
spec="RealRandomWalkOperator" parameter="@seir.a.t:{0}" \
useGaussian="true" weight="1.0" windowSize="5.0"/>'.format(args.out_name)
        template = template.replace(R0Change_operator, '')
        template = template.replace('<!-- R0_CHANGE -->', str(2050))    # sure hope no one is using this code in the year 2050
    else:
        template = template.replace('<!-- R0_CHANGE -->', str(args.R0Change))
    with open(args.base_path+'.xml', 'w') as outfile:
        outfile.write(template)


if __name__ == '__main__':
    generate_xml(args)
