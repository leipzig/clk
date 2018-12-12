import click
import os

@click.command()
@click.argument('manifest')
@click.argument('pair')
@click.argument('cur_dir')
def main(manifest, pair, cur_dir):
    in_bam_str = ''
    in_bam_list = []
    with open(manifest) as in_list:
        line = in_list.readline()
        line = line[:line.index('#')] if '#' in line else line
        n_samp_short = int(line.split()[0])
        for i in range(0, n_samp_short):
            in_bam = []
            line = in_list.readline()
            line = line[:line.index('#')] if '#' in line else line
            n_rep = int(line.split()[0])
            for j in range(0, n_rep):
                line = in_list.readline()
                line = line[:line.index('#')] if '#' in line else line
                in_bam.append(os.path.abspath(os.path.join(cur_dir,line.split()[0])))
            in_bam_list += [in_bam]
            in_bam_str += ','.join(in_bam) if in_bam_str == '' else ':' + ','.join(in_bam)
    print(','.join(in_bam_list[int(pair)-1]))

if __name__ == "__main__":
    main()