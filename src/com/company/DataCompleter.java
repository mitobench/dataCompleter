package com.company;

import com.company.calculations.HaplogrepCaller;
import com.company.calculations.LocationCompleter;
import com.company.calculations.Statistics;
import com.company.io.FastaReader;
import com.company.io.HSDParser;
import com.company.io.MetaInfoReader;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class DataCompleter {

    private String outfile;

    public void run(String data_template_filepath, String data_fasta_filepath, String outfolder) throws IOException {

        String[] fileName = data_fasta_filepath.replaceFirst("[.][^.]+$", "").split("/");
        String fileNameWithoutExt = fileName[fileName.length-1];

        outfile = outfolder + java.time.LocalDateTime.now() + "_" + fileNameWithoutExt + "_" + "_completed.csv";
        BufferedWriter data_meta_file_updated = new BufferedWriter(new FileWriter(outfile));

        MetaInfoReader metaInfoReader = new MetaInfoReader(data_template_filepath);
        metaInfoReader.read();
        String[] types = metaInfoReader.getTypes_list();
        String[] header = metaInfoReader.getHeader_list();

        FastaReader fastaReader = new FastaReader(data_fasta_filepath);
        fastaReader.parseFasta();

        boolean isheaderWritter = false;

        // calculate haplogroups
        HaplogrepCaller haplogrepCaller = new HaplogrepCaller();
        haplogrepCaller.call(data_fasta_filepath);

        HSDParser hsdParser = new HSDParser();
        HashMap<String, ArrayList<String>> entryList = null;

        try {
            hsdParser.parseFile("haplogroups.hsd");
            entryList = hsdParser.getEntryList();
        } catch (Exception e) {
            System.out.println("'haplogroups.hsd' could not be read.");
            System.exit(0);
        }

        haplogrepCaller.deleteTmpFiles();

        // complete geographic locations based on already given information
        LocationCompleter locationCompleter = new LocationCompleter();
        locationCompleter.setHeader(header);
        locationCompleter.setIndexes();

        Statistics calculator = new Statistics();

        ArrayList<String> entries = metaInfoReader.getEntry_list();

        for (String entry : entries) {

            String[] meta_info = entry.split(",", types.length);
            String accessionID_with_version = meta_info[metaInfoReader.getIndexOfArrtibute("##accession_id")].replace("\"","");
            String accessionID = accessionID_with_version.split("\\.")[0].trim();
            String sequence = fastaReader.getSequenceMap().get(accessionID);

            double percentageOfN = calculator.calculatePercentageOfN(sequence);

            // determine user alias
            String user_alias = meta_info[metaInfoReader.getIndexOfArrtibute("user_firstname")].trim() + "" + meta_info[metaInfoReader.getIndexOfArrtibute("user_surname")].trim();
//                String user_alias = meta_info[metaInfoReader.getUserFirstNameIndex()].trim() + "" + meta_info[metaInfoReader.getUserSurnameIndex()].trim();

            // fill HaploGrep2 results
            String haplogroup="NULL";
            String haplotype="NULL";
            String quality="NULL";

            if(entryList == null){
                haplogroup = "NULL";
                haplotype = "NULL";
                quality = "NULL";
            } else {
                try {
                    haplogroup = entryList.get(accessionID).get(0).replace("'", "");
                    haplotype = entryList.get(accessionID).get(3);
                    quality = entryList.get(accessionID).get(1);
                } catch (Exception e) {
                    System.out.println("Sequence with accession id "+ accessionID + " not contained in Haplogrep2 result file");
                }
            }

            // complete geographic information
            locationCompleter.setEntry(entry.split(",", types.length));

            String[] entry_completed = locationCompleter.getCompletedInformation();
            meta_info = entry_completed;
            String meta_info_parsed = "";

            for (int i = 0; i < types.length; i++) {

                String type = types[i].replace("#", "").trim();
                String info;
                if( meta_info[i] == null){
                    info = "NULL";
                } else {
                    info = meta_info[i].replace("\"", "").trim();
                }

                if(info == null) {
                    meta_info_parsed += "NULL,";
                } else if (info.equals("NULL")) {
                    meta_info_parsed += "NULL,";
                }else if (info.equals("")) {
                    meta_info_parsed += "NULL,";
                } else if (info.contains("'")) {
                    String hg_tmp = info.replace("'", "");
                    meta_info_parsed += "'" + hg_tmp + "',";
                } else if (type.equals("String")) {
                    meta_info_parsed += "'" + info + "',";
                } else {
                    meta_info_parsed += info + ",";
                }
            }
            if (meta_info_parsed.endsWith(",")) {
                meta_info_parsed = meta_info_parsed.substring(0, meta_info_parsed.length() - 1);
            }


            // write new header
            if (!isheaderWritter) {
                metaInfoReader.addToHeader(",percentage_N,user_alias,haplogroup_current_versions,macro_haplogroup,haplotype_current_versions,quality_haplotype_current_version, mt_sequence");
                metaInfoReader.addTotypes(",real,String,String,String,String,int,String");
                data_meta_file_updated.write(metaInfoReader.getHeader());
                data_meta_file_updated.newLine();
                data_meta_file_updated.write(metaInfoReader.getTypes());
                data_meta_file_updated.newLine();
                //System.out.println(metaInfoReader.getHeader());
                //System.out.println(metaInfoReader.getTypes());
                isheaderWritter = true;
            }


            // write new meta data entry
            String values = meta_info_parsed + "," + percentageOfN + ",'" + user_alias + "','" + haplogroup + "','" +
                    haplogroup.substring(0,2) + "," + haplotype + "'," + quality + ",'" + fastaReader.getSequenceMap().get(accessionID) + "'";

            values = values.replace("'NULL'", "NULL");
            values = values.replace("NULL", "");
            values = values.replace("'\"", "'");
            values = values.replace("\"'", "'");
            values = values.replace("\"", "");
            values = values.replace("'", "");

            data_meta_file_updated.write(values);
            data_meta_file_updated.newLine();
            //System.out.println(values);
            //System.out.println("");
        }

        data_meta_file_updated.close();
    }

    public String getOutfile() {
        return outfile;
    }

}
