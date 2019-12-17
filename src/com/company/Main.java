package com.company;

import java.io.IOException;

public class Main {

    public static void main(String[] args) throws IOException {

        String help = "";

        if(args.length < 2){
            System.out.println(help);
            System.exit(0);

        } else {

            String data_template_filepath = args[0];    // mandatory
            String data_fasta_filepath = args[1];       // mandatory
            String outfolder="";
            if(args.length == 3)
                outfolder = args[2];                 // optional

            DataCompleter runner = new DataCompleter();
            runner.run(data_template_filepath, data_fasta_filepath, outfolder);

        }
    }
}
