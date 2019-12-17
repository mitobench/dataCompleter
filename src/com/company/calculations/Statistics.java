package com.company.calculations;

public class Statistics {

    public double calculatePercentageOfN(String sequence) {

        int count_n = 0;
        for(char c : sequence.toCharArray()){
            if(c == 'N'){
                count_n++;
            }
        }

        return count_n/(double)sequence.length();

    }

    public double calculateCompleteness(String sequence) {
        return sequence.length() / 16569.0;
    }
}
