#! /usr/bin/python

# produce one ptt / record present in the genbank file

from Bio import SeqIO
from optparse import OptionParser
import re
from Bio import SeqUtils

def circos_fasta_draft(fasta_file_name):
    from Bio import SeqIO
    contigs = []
    start = 1
    handle = open(fasta_file_name)
    for seq_record in SeqIO.parse(handle, "fasta"):
        stop = start + len(seq_record.seq)
        contigs.append([seq_record.name, start, stop])
        start = stop+1
    return contigs


def circos_fasta_draft_misc_features(record):
    gap_locations = []
    for feature in record.features:
        if feature.type == "assembly_gap":
            gap_locations.append(feature.location)
    if len(gap_locations) == 0:
        return None
    contigs = []
    for i in range(0, len(gap_locations)):
        if i == 0:
            contigs.append([record.name + "_%s" % 1, 0, int(gap_locations[i].start)])
        else:
            contigs.append([record.name + "_%s" % (i+1), int(gap_locations[i-1].end), int(gap_locations[i].start)])
    contigs.append([record.name + "_%s" % (i+2), int(gap_locations[-1].end), int(len(record.seq))])
    return contigs


def find_index(pattern, seq):
  """Return first item in sequence where f(item) == True."""
  for item in seq:
    if re.match(pattern,item):
      return seq.index(item)


def print_circos_record_file(record_list, out_name = "circos_contigs.txt", draft_data = False, draft_coordinates=False):
    i = 1
    f = open(out_name, "w")
    x = 0

    for record in record_list:

      try:
          for contig in draft_data[x]:
              if contig[2] < contig[1]:
                  #print 'impossible contig limits!', contig
                  continue
              #print contig
              if contig[1] == contig[2]:
                  continue
              if draft_coordinates:
                  if i%2 == 0:
                      f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (contig[0], contig[0], 0, contig[2]-contig[1], 3))
                  else:
                      ##print "chr -", record.contig, record.contig, "0",len(record.seq), "spectral-5-div-%s" % (4)
                      f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (contig[0], contig[0], 0, contig[2]-contig[1], 4))
                  i+=1
              else:
                  if i%2 == 0:
                      f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (contig[0], contig[0], contig[1], contig[2], 3))
                  else:
                      ##print "chr -", record.contig, record.contig, "0",len(record.seq), "spectral-5-div-%s" % (4)
                      f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (contig[0], contig[0], contig[1], contig[2], 4))
                  i+=1

          x+=1
      except:
            if i%2 == 0:
                f.write("chr - %s %s 0 %s spectral-5-div-%s\n" % (record.id.split(".")[0], record.id.split(".")[0], len(record.seq), 3))
            else:
                ##print "chr -", record.contig, record.contig, "0",len(record.seq), "spectral-5-div-%s" % (4)
                f.write("chr - %s %s 0 %s spectral-5-div-%s\n" % (record.id.split(".")[0], record.id.split(".")[0], len(record.seq), 4))
            i += 1
            x += 1
    f.close()


def print_circos_labels_file(record_list, locus2label, out_name = "circos_labels.txt",
                           draft_data=[None], draft_coordinates=False):

    import numpy

    f = open(out_name, "w")

    for y, record in enumerate(record_list):
        for feature in record.features:
            if feature.type == "CDS":

                    try:
                        for i in draft_data[y]:
                            # determine to which contig the feature belong

                            if feature.location.start >= i[1] and feature.location.end <= i[2]:
                                if draft_coordinates:
                                    contig = i[0]
                                    start = feature.location.start - i[1]
                                    end = feature.location.end - i[1]
                                else:
                                    contig = i[0]
                                    start = feature.location.start
                                    end = feature.location.end
                    except IndexError:
                        contig = record.id.split(".")[0] # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    except TypeError:
                        contig = record.id.split(".")[0] # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end

                    # absurd features
                    if numpy.abs(feature.location.start-feature.location.end) > 50000:
                        continue

                    if 'pseudo' in feature.qualifiers:
                        continue

                    if feature.qualifiers['orthogroup'][0] in locus2label or feature.qualifiers['locus_tag'][0] in locus2label:

                            f.write('%s %s %s %s color=spectral-5-div-2,z=1\n' % (contig,
                                                                               start,
                                                                               end,
                                                                               locus2label[feature.qualifiers['locus_tag'][0]]))

    f.close()


def print_circos_gene_file(record_list, feature_type="CDS", strand ="1",
                           out_name = "circos_genes_plus.txt",
                           locus_highlight=[],
                           draft_data=[None],
                           group_id2orthologs_presence=False,
                           query_taxon_id=False,
                           color_missing=True,
                           draft_coordinates=False,
                           locus_highlight2=[]):

    #print "highlight:", locus_highlight

    import numpy

    if strand == "1":
        f = open(out_name, "w")
    if strand == "-1":
        f = open(out_name, "w")

    #print 'n records:', len(record_list)

    for y, record in enumerate(record_list):
        for feature in record.features:
            if feature.type == feature_type:
                if str(feature.strand) == strand:
                    try:
                        for i in draft_data[y]:
                            # determine to which contig the feature belong

                            if feature.location.start >= i[1] and feature.location.end <= i[2]:
                                if draft_coordinates:
                                    contig = i[0]
                                    start = feature.location.start - i[1]
                                    end = feature.location.end - i[1]
                                else:
                                    contig = i[0]
                                    start = feature.location.start
                                    end = feature.location.end
                    except IndexError:
                        contig = record.id.split(".")[0] # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    except TypeError:
                        contig = record.id.split(".")[0] # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end

                    '''
                    # in case the second record is not fragmented in several contigs (no draft data)
                    except IndexError:
                        #print "no draft for", record.id
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    except TypeError:
                        #print "no draft for", record.id
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    '''
                    #print "longueur", len(feature.location), type(len(feature.location)), feature.location.start, feature.location.end


                    if numpy.abs(feature.location.start-feature.location.end) > 50000:
                        continue

                    if 'pseudo' in feature.qualifiers:

                        continue
                    try:
                        a = feature.qualifiers['orthogroup']
                        orthology_tag = True
                    except:
                        orthology_tag = False

                    if orthology_tag:
                        if feature.qualifiers['orthogroup'][0] in locus_highlight or feature.qualifiers['locus_tag'][0] in locus_highlight:

                            #print "COLORS!!!!!!!!!!!!!!!!!!!"
                            try:
                                '''
                                f.write('%s %s %s fill_color=spectral-5-div-4,id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                                         start,
                                                                                                         end,
                                                                                                         feature.qualifiers['locus_tag'][0],
                                                                                                         feature.qualifiers['gene'][0],
                                                                                                         re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                                         feature.qualifiers['product'][0])))
                                '''
                                f.write('%s %s %s fill_color=piyg-5-div-1,id=%s,z=1\n' % (contig,
                                                                        start,
                                                                        end, feature.qualifiers['locus_tag'][0]))


                            except:
                                '''
                                f.write('%s %s %s fill_color=spectral-5-div-4,id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                                         start,
                                                                                                         end,
                                                                                                         feature.qualifiers['locus_tag'][0],
                                                                                                         "-",
                                                                                                         re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                                         feature.qualifiers['product'][0])))
                                '''
                                f.write('%s %s %s fill_color=spectral-5-div-4,id=%s,z=1\n' % (contig,
                                                                        start,
                                                                        end,
                                                                        feature.qualifiers['locus_tag'][0]))
                        elif feature.qualifiers['orthogroup'][0] in locus_highlight2 or feature.qualifiers['locus_tag'][0] in locus_highlight2:
                            try:
                                f.write('%s %s %s fill_color=paired-11-qual-10,id=%s,z=1\n' % (contig,
                                                                        start,
                                                                        end, feature.qualifiers['locus_tag'][0]))


                            except:
                                f.write('%s %s %s fill_color=spectral-5-div-4,id=%s,z=1\n' % (contig,
                                                                        start,
                                                                        end,
                                                                        feature.qualifiers['locus_tag'][0]))


                        else:

                            f.write('%s %s %s id=%s\n' % (contig,
                                                  start,
                                                  end, feature.qualifiers['locus_tag'][0]))

                    elif feature.qualifiers['locus_tag'][0] in locus_highlight:
                        try:

                            f.write('%s %s %s fill_color=piyg-5-div-1,id=%s,z=1\n' % (contig,
                                                                    start,
                                                                    end, feature.qualifiers['locus_tag'][0]))


                        except:
                            f.write('%s %s %s fill_color=spectral-5-div-4,id=%s,z=1\n' % (contig,
                                                                    start,
                                                                    end,
                                                                    feature.qualifiers['locus_tag'][0]))

                        else:

                            f.write('%s %s %s id=%s\n' % (contig,
                                                  start,
                                                  end, feature.qualifiers['locus_tag'][0]))
                    # not in locus_highlight
                    else:
                        if group_id2orthologs_presence and query_taxon_id and color_missing:
                            try:
                                pseudo = feature.qualifiers['pseudogene']
                                #print 'pseudogene, continue'
                                continue
                            except:
                                pass
                            if group_id2orthologs_presence[feature.qualifiers['orthogroup'][0]][int(query_taxon_id)] == 0:
                                f.write('%s %s %s fill_color=orange\n' % (contig,
                                                                        start,
                                                                        end))
                            else:
                                f.write('%s %s %s id=%s\n' % (contig,
                                                      start,
                                                      end, feature.qualifiers['locus_tag'][0]))
                        else:

                            f.write('%s %s %s id=%s\n' % (contig,
                                                  start,
                                                  end, feature.qualifiers['locus_tag'][0]))

                        '''
                            try:
                                f.write('%s %s %s fill_color=orange, id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                                             start,
                                                                                                             end,
                                                                                                             feature.qualifiers['locus_tag'][0],
                                                                                                             feature.qualifiers['gene'][0],
                                                                                                             re.sub("[ |\-|(|)|\]|\[|\.|,]+", "_",
                                                                                                                    feature.qualifiers['product'][0])))
                            except:
                                f.write('%s %s %s fill_color=orange, id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                                             start,
                                                                                                             end,
                                                                                                             feature.qualifiers['locus_tag'][0],
                                                                                                             "-",
                                                                                                             re.sub("[ |\-|(|)|\]|\[|\.|,]+", "_",
                                                                                                                    feature.qualifiers['product'][0])))

                        else:

                            try:
                                f.write('%s %s %s id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                     start,
                                                                                     end,
                                                                                     feature.qualifiers['locus_tag'][0],
                                                                                     feature.qualifiers['gene'][0],
                                                                                     re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                            feature.qualifiers['product'][0])))
                            except:
                                f.write('%s %s %s id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                     start,
                                                                                     end,
                                                                                     feature.qualifiers['locus_tag'][0],
                                                                                     '-',
                                                                                     re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                            feature.qualifiers['product'][0])))
                    else:

                        try:
                            f.write('%s %s %s id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                 start,
                                                                                 end,
                                                                                 feature.qualifiers['locus_tag'][0],
                                                                                 feature.qualifiers['gene'][0],
                                                                                 re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                        feature.qualifiers['product'][0])))
                        except:
                            try:
                                f.write('%s %s %s id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                 start,
                                                                                 end,
                                                                                 feature.qualifiers['locus_tag'][0],
                                                                                 '-',
                                                                                 re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                        feature.qualifiers['product'][0])))
                            except:
                                f.write('%s %s %s id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                 start,
                                                                                 end,
                                                                                 feature.qualifiers['locus_tag'][0],
                                                                                 '-',
                                                                                 re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                        "-")))
                        '''

            if feature.type == 'rRNA':
                #print 'rrna--------------------'
                if str(feature.strand) == strand:
                    try:
                        for i in draft_data[y]:
                            # determine to which contig the feature belong

                            if feature.location.start >= i[1] and feature.location.end <= i[2]:
                                if draft_coordinates:
                                    contig = i[0]
                                    start = feature.location.start - i[1]
                                    end = feature.location.end - i[1]
                                else:
                                    contig = i[0]
                                    start = feature.location.start
                                    end = feature.location.end
                    except IndexError:
                        contig = record.id.split(".")[0] # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    except TypeError:
                        contig = record.id.split(".")[0] # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end

                    #print 'rrna position:', contig, start, end, feature.qualifiers
                    f.write('%s %s %s fill_color=pblue\n ' % (contig,
                                                            start,
                                                            end))


            if feature.type == 'tRNA':
                if str(feature.strand) == strand:
                    try:
                        for i in draft_data[y]:
                            # determine to which contig the feature belong

                            if feature.location.start >= i[1] and feature.location.end <= i[2]:
                                if draft_coordinates:
                                    contig = i[0]
                                    start = feature.location.start - i[1]
                                    end = feature.location.end - i[1]
                                else:
                                    contig = i[0]
                                    start = feature.location.start
                                    end = feature.location.end
                    except IndexError:
                        contig = record.id.split(".")[0] # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    except TypeError:
                        contig = record.id.split(".")[0] # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end

                    f.write('%s %s %s fill_color=pred\n ' % (contig,
                                                            start,
                                                            end))


    f.close()


def print_circos_GC_file(record_list, 
                         feature_type="CDS", 
                         out_directory="",
                         outname=False):

    import GC
    import os

    if not outname:
        out_var_file = os.path.join(out_directory, 'circos_GC_var_%s.txt' % record_list[0].name)
        out_skew_file = os.path.join(out_directory, 'circos_GC_skew_%s.txt' % record_list[0].name)
    else:
        out_var_file, out_skew_file = outname
        
    f = open(out_var_file, 'w')
    g = open(out_skew_file, 'w')

    out_var = ''
    out_skew = ''
    for record in record_list:
        # this function handle scaffolds (split sequence when encountering NNNNN regions)
        out_var += GC.circos_gc_var(record, windows=1000)
        out_skew += GC.circos_gc_skew(record, windows=1000)
    f.write(out_var)
    g.write(out_skew)

    return (out_var_file, out_skew_file)



class Circos_config:
  def __init__(self, caryotype_file,
               chr_spacing_list=[],
               show_ticks="yes",
               show_tick_labels="yes",
               ideogram_spacing=0,
               radius=0.75,
               label_radius=0.05,
               show_ideogram_labels="no",
               color_files=""):
    import re

    self.plots = ""
    self.links = ""
    self.highlights = ""

    self.template_caryotype= "karyotype = %s\n" \
                             " chromosomes_units           = 10000\n" \
                             " chromosomes_display_default = yes\n" % caryotype_file


    self.template_ideograms = "<ideogram>\n" \
                              " <spacing>\n" \
                              " default            = %su\n" \
                              " %s" \
                              " </spacing>\n" \
                              " \n" \
                              " # thickness and color of ideograms\n" \
                              " thickness          = 12p\n" \
                              " stroke_thickness   = 1\n" \
                              " stroke_color       = black\n" \
                              " \n" \
                              " # the default chromosome color is set here and any value\n" \
                              " # defined in the karyotype file overrides it\n" \
                              " fill               = yes\n" \
                              " fill_color         = black\n" \
                              " \n" \
                              " # fractional radius position of chromosome ideogram within image\n" \
                              " radius             = %sr\n" \
                              " show_label         = %s\n" \
                              " label_font         = default\n" \
                              " label_radius       = dims(ideogram,radius) + %sr\n" \
                              " label_size         = 30\n" \
                              " label_parallel     = no\n" \
                              " \n" \
                              " # show_bands determines whether the outline of cytogenetic bands\n" \
                              " # will be seen\n" \
                              " show_bands         = yes\n" \
                              " band_stroke_thickness = 1\n" \
                              " " \
                              " # in order to fill the bands with the color defined in the karyotype\n" \
                              " # file you must set fill_bands\n" \
                              " fill_bands         = yes\n" \
                              " band_transparency  = 1\n" \
                              " \n" \
                              " </ideogram>\n" % (ideogram_spacing,
                                                  self.add_spacing(chr_spacing_list),
                                                  radius,
                                                  show_ideogram_labels,
                                                  label_radius)

    self.end_ticks = "  <tick>\n" \
                      "  multiplier   = 1\n" \
                      "  position = end\n" \
                      " show_ticks         = yes\n" \
                      "  size              = 4p\n" \
                      "  show_label        = no\n" \
                      "  label_size        = 0p\n" \
                      "  format    = %s bp\n" \
                      "  </tick>\n" \

    self.big_ticks = " <tick>\n" \
                          " show_ticks         = yes\n" \
                          " skip_first_label = no\n" \
                          " multiplier   = 10/1u\n" \
                          " spacing           = 10u\n" \
                          " size              = 15p\n" \
                          " show_label        = yes\n" \
                          " label_size        = 25p\n" \
                          " format            = %s kb\n" \
                          " thickness         = 2p\n" \
                          " </tick>\n" \


    self.template_ticks = "show_ticks         = %s\n" \
                          " show_tick_labels   = %s\n" \
                          " \n" \
                          " <ticks>\n" \
                          " tick_label_font    = condensed\n" \
                          " radius             = dims(ideogram,radius_outer)\n" \
                          " label_offset       = 8p\n" \
                          " label_size         = 4p\n" \
                          " color              = black\n" \
                          " thickness          = 2p\n" \
                          " \n" \
                          " \n" \
                          " <tick>\n" \
                          " show_ticks         = yes\n" \
                          " skip_first_label = no\n" \
                          " multiplier   = 10/1u\n" \
                          " spacing           = 1u\n" \
                          " size              = 5p\n" \
                          " show_label        = no\n" \
                          " label_size        = 5p\n" \
                          " format            = %s\n" \
                          " </tick>\n" \
                          " %s\n" \
                          " \n" \
                          " \n" \
                          " </ticks>\n" % (show_ticks, show_tick_labels, "%%.1d", self.big_ticks % "%%.1d")





    self.template_rules = "<rules>\n" \
                          " %s\n" \
                          "</rules>\n"

    self.template_backgrounds = "<backgrounds>\n" \
                          " %s\n" \
                          "</backgrounds>\n"

    # #" <<include colors.rn.conf>>\n" \
    self.settings ="<colors>\n" \
                   " %s\n" \
                   " #<<include brewer.all.conf>>\n" \
                   " </colors>\n" \
                   " <image>\n"\
                   " image_map_use      = yes\n" \
                   " image_map_overlay  = no\n" \
                   " image_map_overlay_stroke_color     = red\n" \
                   " <<include image.conf>>\n" \
                   " </image>\n" \
                   " #includes  etc/colors.conf\n" \
                   " #          etc/fonts.conf\n" \
                   " #          etc/patterns.conf\n" \
                   " <<include colors_fonts_patterns.conf>>\n" \
                   " # system and debug settings\n" \
                   " <<include housekeeping.conf>>\n" \
                   " anti_aliasing*     = no\n"

    self.complete_file = self.template_caryotype + self.template_ideograms + self.template_ticks + "%s %s %s" + self.settings % (color_files)

  def _template_spacing(self, chr1, chr2):
    template = '<pairwise %s %s>\n' \
               ' spacing = 2u\n' \
               '</pairwise>\n' % (chr1, chr2)
    return template

  def _template_plot(self, file, type="line", r0=1,
                     r1=1.05, color=False, fill_color="red", thickness = "0.8p", z = 1, rules ="", backgrounds="", url="", min=False, max=False):

        #print 'template color------------------', color
        template1 = "<plot>\n" \
               "type		    = %s\n" \
               " url                = %s[id]\n" \
               " r0                 = %s\n" \
               " r1                 = %s\n" \
               " fill_color         = %s\n" \
               " thickness          = %s\n" \
               " file               = %s\n" \
               " z                  = %s\n" % (type, url, r0, r1, fill_color, thickness, file, z)
        #print '--------------- COLOR--------------'
        #print color
        if color:
            template1+= " color          = %s\n" % color
        if min:
            template1+= " min          = %s\n" % min
        if max:
            max+= " max          = %s\n" % max
        template_rules = " %s\n" \
               " %s\n" \
               " </plot>\n" % (rules, backgrounds)

        return template1 + template_rules

  def template_rule(self, condition, fill_color):
    template = "<rule>\n" \
               " condition          = %s\n" \
               " fill_color         = %s\n" \
               " </rule>\n" % (condition, fill_color)
    return template

  def template_background(self, color):
    template = "<background>\n" \
               " color         = %s\n" \
               " </background>\n" % (color)
    return template



  def _template_link(self, link_file, color="black_a5", thickness=1):
    template = "<link>\n" \
               "ribbon = yes\n" \
               "file          = %s\n" \
               "color         = %s\n" \
               "radius        = 0.99r\n" \
               "bezier_radius = 0.1r\n" \
               "thickness     = %s\n" \
               " </link>\n" % (link_file, color, thickness)
    return template


  def _template_highlight(self, file, fill_color="grey_a1", r1="1.55r", r0="1.50r", url="/chlamdb/locusx/chlamydia_03_15/"):
    template ="<highlight>\n" \
              " fill_color = %s\n" \
              " file       = %s\n" \
              " r1         = %s\n" \
              " r0         = %s\n" \
              " url = %s[id]\n" \
              " </highlight>\n" % (fill_color, file, r1, r0, url)

    return template


  def add_plot(self, file, type="line", r0=1, r1=1.05,
               fill_color="grey_a1", thickness = "2p", z = 1, rules ="", backgrounds="", url="", min=False, max=False, color=False):
    plot = self._template_plot(file, type, r0, r1, color, fill_color, thickness, z, rules, backgrounds, url, min=min, max=max)
    if len(re.findall("</plots>", self.plots))>0:
      # remove end balise
      self.plots = re.sub("</plots>", "", self.plots)
      # add new plot and end balise
      self.plots = self.plots + plot + "</plots>\n"
    else:
      self.plots = "<plots>\n" + plot + "</plots>\n"



  def add_link(self, link_file, color="black_a5", thickness=1):
    link = self._template_link(link_file, color=color, thickness=thickness)
    if len(re.findall("</links>", self.plots))>0:
      # remove end balise
      self.links = re.sub("</links>", "", self.plots)
      # add new plot and end balise
      self.links = self.links + link + "</links>\n"
    else:
      self.links = "<links>\n" + link + "</links>\n"


  def add_highlight(self, file, fill_color="grey_a1", r1="1.55r", r0="1.50r", href=""):
    highlight = self._template_highlight(file, fill_color, r1, r0, href)
    if len(re.findall("</highlights>", self.highlights))>0:
      # remove end balise
      self.highlights = re.sub("</highlights>", "", self.highlights)
      # add new plot and end balise
      self.highlights = self.highlights + highlight + "</highlights>\n"
    else:
      self.highlights = "<highlights>\n" + highlight + "</highlights>\n"

  def add_spacing(self, chr_pair_list):
      all_spacing = ''
      if len(chr_pair_list) == 0:
          return ''
      else:
          for pair in chr_pair_list:
              all_spacing += self._template_spacing(pair[0], pair[1])
          return all_spacing

  def get_file(self):
    return self.complete_file % (self.plots, self.highlights, self.links)


def get_circos_GC_config_files(biodatabase_name, accession_list):
    from chlamdb.plots import GC
    '''
    accessions: in case of several chromosomes or combinations of chromosomes

    '''

    server, db = manipulate_biosqldb.load_db(biodatabase_name)

    final_gc_var = ''
    final_gc_skew = ''

    for accession in accession_list:
        record = db.lookup(accession=accession)

        #print record
        mean_GC = GC.GC(record.seq)
        #print "mean GC", mean_GC
        circos_gc_var = GC.circos_gc_var(record)
        circos_gc_skew = GC.circos_gc_skew(record)
        final_gc_var += circos_gc_var
        final_gc_skew += circos_gc_skew

        #import pylab
        #import numpy as np
        #cumulated_skew = np.cumsum(values)
        #pylab.plot(cumulated_skew)
        #pylab.plot(values[0:-1])
        #pylab.show()

    return (final_gc_var, final_gc_skew)


