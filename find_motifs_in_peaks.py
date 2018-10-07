import sys

'''
Program is to correlate the identified motifs from RSAT to which peaks 
the motif resides in and output that to a table that can be sorted and viewed.

Author: Debbie Thurtle-Schmidt
Date: 10/19/2017

'''

def read_peaks(peak_file):

    """
    reads in a peak file in bed format (ex. from MACS2)

    Input: peak file name/path

    Return: a list of each line, with each line as a list

    """

    with open(peak_file, 'rU') as fh1:
        peak_list = []
        for line in fh1:
            peak_list.append(line.strip().split('\t'))

    return peak_list


def read_motifs(motif_file):
    '''

    reads in a motif file from the tab output of RSAT

    Args:
    tab output from RSAT file path

    Returns:
    a list of only the lines that have motif coordinates, with each
    motif coordinate as a list.

    '''
    with open(motif_file, 'rU') as fh2:
        motif_list = []
        for line in fh2:
            if line[0].startswith(';') or line[0].startswith('#'):
                pass
            else:
                motif_list.append(line.strip().split('\t'))

    return motif_list


def motifs_in_peaks(peak_list, motif_list, motif_table_layout):
    '''

    identifies which peaks have which motifs. For each peak, it will be recorded
    if it has that motif.

    Input: peak_list - a list of peaks
           motif_list - a list of motifs
           motif_table - empty list of the correct length for the number of columns
           in the table
           motif_table_layout - a dictionary where the value is the name of the
           motif and the key is the index that motif corresponds to when recording
           if the peak has that motif

    Return: a list of each peak with each motif it has
    '''
    peak_motif_table = []
    for peak in peak_list:
        peak_location = peak[0] + ':' + peak[1] + '-' + peak[2] #make string that looks like motif_list[0]
        peak_row = ['0'] * 12 #using magic constant for now
        peak_row[0] = peak_location
        for motif in motif_list:
            #append the weight of the motif to the correct index in the list
            if motif[0] == peak_location:
                print motif[0], motif[2], peak_location
                index = motif_table_layout[motif[2]]
                #if peak_row[index] == 'NA':
                peak_row[index] = motif[7]
                #else:
                    #peak_row[index] = peak_row[index] + ',' + motif[7]
            else:
                pass
        print peak_row
        peak_motif_table.append(peak_row)

    return peak_motif_table

def correct_table(motif_table):
    '''

    Any place there is a blank spot (that motif doesn't exist) in that list,
    change it to a NA.

    Input: a list of lists creating the peaks with each motif found in that peak

    Return: a list of each peak with each motif. If a motif is not found in that
    peak put NA

    '''

def main():
    peak_list = read_peaks(sys.argv[1])
    motif_list = read_motifs(sys.argv[2])
    #make a dictionary to input into motifs_in_peaks
    '''
    motif_dictionary = {'oligos_6nt_mkv3_m2' : 1, 'oligos_7nt_mkv4_m1' : 1,
                         'dyads_m6' : 1, 'dyads_m1' : 1, 'dyads_m10' : 1, 'dyads_m9' : 1,
                         'oligos_7nt_mkv4_m6' : 2, 'oligos_7nt_mkv4_m4' : 2, 'oligos_7nt_mkv4_m2' : 2,
                         'oligos_6nt_mkv3_m1' : 2, 'oligos_6nt_mkv3_m7' : 2, 'local_words_6nt' : 1,
                         'positions_7nt_m6' : 1, 'oligos_7nt_mkv4_m3' : 3, 'oligos_7nt_mkv4_m5' : 3,
                         'oligos_7nt_mkv4_m9' : 4, 'dyads_m3' : 4, 'dyads_m2' : 4,
                         'dyads_m4' : 4, 'dyads_m5' : 4, 'oligos_6nt_mkv3_m8' : 4, 'dyads_m8' : 5,
                         'dyads_m7' : 5, 'oligos_6nt_mkv3_m3' : 5, 'oligos_7nt_mkv4_m10' : 6,
                         'oligos_6nt_mkv3_m9' : 6, 'oligos_6nt_mkv3_m10' : 7, 'oligos_6nt_mkv3_m5' : 7,
                         'positions_7nt_m1' : 8, 'positions_7nt_m5' : 8, 'oligos_7nt_mkv4_m7' : 11,
                         'positions_7nt_m2' : 11, 'positions_7nt_m3' : 11, 'positions_7nt_m4' : 11,
                         'oligos_6nt_mkv3_m4' : 11, 'oligos_7nt_mkv4_m8' : 9,
                         'positions_7nt_m7' : 11, 'oligos_6nt_mkv3_m6' : 10}
                         '''
    motif_dictionary = {'dyads_m9' : 1, 'local_words_6nt' : 1, 'positions_6nt' : 1,
                        'dyads_m10' : 1, 'oligos_7nt_mkv4_m1' : 1, 'oligos_6nt_mkv3_m1' : 1,
                        'dyads_m1' : 1, 'dyads_m7' : 1, 'oligos_7nt_mkv4_m6' : 2,
                        'oligos_7nt_mkv4_m4' : 2, 'oligos_7nt_mkv4_m6' : 3, 'oligos_7nt_mkv4_m7' : 3,
                        'oligos_7nt_mkv4_m2' : 3, 'oligos_6nt_mkv3_m2' : 3,
                        'dyads_m2' : 4, 'dyads_m4' : 4, 'dyads_m3' : 4, 'dyads_m5' : 4,
                        'dyads_m6' : 4, 'dyads_m8' : 4, 'oligos_6nt_mkv3_m3' : 4,
                        'oligos_7nt_mkv4_m10' : 5, 'oligos_6nt_mkv3_m9' : 5, k}
    peak_motif_table = motifs_in_peaks(peak_list, motif_list, motif_dictionary)
    motif_header = ['peak_name', 'cluster_1', 'cluster_2' , 'cluster_3', 'cluster_4',
                    'cluster_5', 'nhr-48', 'Y67B8A.3/fkh-6',
                    'irx-1/daf-16', 'nhr-23', 'klf-2', 'rest']
    peak_motif_table_file = []
    for entry in peak_motif_table:
        peak_motif_table_file.append('\t'.join(entry) + '\n')
    with open(sys.argv[1] + '.motifs_zeros_fewermotifs', 'w') as fh_write:
        fh_write.write('\t'.join(motif_header) + '\n')
        fh_write.writelines(peak_motif_table_file)


if __name__ == "__main__":
    main()
