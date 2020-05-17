#!/usr/bin/env python3
from Bio import SeqIO
import random

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
    template = template.replace('<!-- EM -->', str(args.EM))
    template = template.replace('<!-- IhM -->', str(args.IhM))
    template = template.replace('<!-- EXOG_PERIOD -->',
                                str(365/(args.infPeriod+args.latPeriod)))

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
    if args.R0ChangeDate == 'NA':
        # todo make this...better?
        R0ChangeDate_operator = '''<operator id="seir.a.operator.t:{0}" \
spec="RealRandomWalkOperator" parameter="@seir.a.t:{0}" \
useGaussian="true" weight="1.0" windowSize="0.5"/>'''.format(args.out_name)
        template = template.replace(R0ChangeDate_operator, '')
        template = template.replace('<!-- R0_CHANGE -->', str(2050))    # sure hope no one is using this code in the year 2050
    else:
        template = template.replace('<!-- R0_CHANGE -->', str(args.R0ChangeDate))
    with open(args.base_path+'/'+args.out_name+'.xml', 'w') as outfile:
        outfile.write(template)


if __name__ == '__main__':
    generate_xml(args)
