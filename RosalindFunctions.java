import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import static java.util.Map.entry;

public class RosalindFunctions {

    final static Map<String, String> CODON_TABLE = Map.ofEntries(
        entry("UUU", "F"), entry("CUU", "L"), entry("AUU", "I"), entry("GUU", "V"),
        entry("UUC", "F"), entry("CUC", "L"), entry("AUC", "I"), entry("GUC", "V"),
        entry("UUA", "L"), entry("CUA", "L"), entry("AUA", "I"), entry("GUA", "V"),
        entry("UUG", "L"), entry("CUG", "L"), entry("AUG", "M"), entry("GUG", "V"),
        entry("UCU", "S"), entry("CCU", "P"), entry("ACU", "T"), entry("GCU", "A"),
        entry("UCC", "S"), entry("CCC", "P"), entry("ACC", "T"), entry("GCC", "A"),
        entry("UCA", "S"), entry("CCA", "P"), entry("ACA", "T"), entry("GCA", "A"),
        entry("UCG", "S"), entry("CCG", "P"), entry("ACG", "T"), entry("GCG", "A"),
        entry("UAU", "Y"), entry("CAU", "H"), entry("AAU", "N"), entry("GAU", "D"),
        entry("UAC", "Y"), entry("CAC", "H"), entry("AAC", "N"), entry("GAC", "D"),
        entry("UAA", ""), entry("CAA", "Q"), entry("AAA", "K"), entry("GAA", "E"),
        entry("UAG", ""), entry("CAG", "Q"), entry("AAG", "K"), entry("GAG", "E"),
        entry("UGU", "C"), entry("CGU", "R"), entry("AGU", "S"), entry("GGU", "G"),
        entry("UGC", "C"), entry("CGC", "R"), entry("AGC", "S"), entry("GGC", "G"),
        entry("UGA", ""), entry("CGA", "R"), entry("AGA", "R"), entry("GGA", "G"),
        entry("UGG", "W"), entry("CGG", "R"), entry("AGG", "R"), entry("GGG", "G")
    );

    public static void countingNucleotides(String filename) {
        try {
            FileReader reader = new FileReader(filename);
            BufferedReader br = new BufferedReader(reader);

            int A = 0;
            int T = 0;
            int C = 0;
            int G = 0;

            String line = br.readLine();

            while (line != null) {
                for (int i = 0; i < line.length(); i++) {
                    switch (line.charAt(i)) {
                        case 'A':
                            A++;
                            break;
                        case 'T':
                            T++;
                            break;
                        case 'C':
                            C++;
                            break;
                        case 'G':
                            G++;
                            break;
                        default:
                            break;
                    }

                }
                line = br.readLine();
            }
            System.out.print(A + " ");
            System.out.print(C + " ");
            System.out.print(G+ " ");
            System.out.print(T + " ");

            br.close();




        } catch (FileNotFoundException e) {
            System.out.println("File not found");
        } catch (IOException e) {
            System.out.println("Error reading file");
        }
    }

    public static void DnaToRna(String filename) {
        try {
            FileReader reader = new FileReader(filename);
            BufferedReader br = new BufferedReader(reader);

            String line = br.readLine();

            String rna = "";

            while (line != null) {
                for (int i = 0; i < line.length(); i++) {
                    if (line.charAt(i) == 'T') {
                        rna = rna + "U";
                    } else {
                        rna = rna + line.charAt(i);
                    }
                }
                line = br.readLine();
            }
            System.out.println(rna);
            br.close();

        } catch (FileNotFoundException e) {
            System.out.println("File not found");
        } catch (IOException e) {
            System.out.println("Error reading file");
        }
    
    }

    public static void complementStrand(String filename) {

        try {
            FileReader reader = new FileReader(filename);
            BufferedReader br = new BufferedReader(reader);
            String line = br.readLine();
            StringBuilder firstStrand = new StringBuilder();
            String complementStrand = "";
            while (line != null) {
                for (int i = 0; i < line.length(); i++) {
                    firstStrand.append(line.charAt(i));
                }
                line = br.readLine();
            }
            String firstStrandString = firstStrand.reverse().toString();
            for (int i = 0; i < firstStrandString.length(); i++) {
                switch (firstStrandString.charAt(i)) {
                    case 'A':
                        complementStrand = complementStrand + "T";
                        break;
                        case 'C':
                            complementStrand = complementStrand + "G";
                            break;
                            case 'G':
                                complementStrand = complementStrand + "C";
                                break;
                                case 'T':
                                    complementStrand = complementStrand + "A";
                                    break;
                                    default:
                                        break;
                }
            }
            System.out.println(complementStrand);
            br.close();


        } catch (FileNotFoundException e) {
            System.out.println("File not found");
        } catch (IOException e) {
            System.out.println("Error reading file");
        }

    }

    public static void maxGCContent(String filename) {

        try {

            FileReader reader = new FileReader(filename);
            BufferedReader br = new BufferedReader(reader);

            HashMap<String, String> sequences = new HashMap<>();
            HashMap<String, Float> gcContent = new HashMap<>();
            String id = "";
            String sequence = "";

            String line = br.readLine();

            while (line != null ) {

                if (line.startsWith(">")) {
                    if (id != "" && sequence != "") id = id.substring(1);
                    sequences.put(id, sequence);
                    id = line;
                    sequence = "";
                } else {
                   sequence += line;
                }


                line = br.readLine();
            }


            for (String key : sequences.keySet()) {
                String seq = sequences.get(key);
                int gcCount = 0;
                for (int i = 0; i < seq.length(); i++) {
                    if (seq.charAt(i) == 'G' || seq.charAt(i) == 'C') {
                        gcCount++;
                    }
                }
                float gcPercentage = ((float) gcCount / seq.length()) * 100;
                gcContent.put(key, gcPercentage);

            }
            
            float currentMax = Integer.MIN_VALUE;
            String maxContentId = "";

            for (String key : gcContent.keySet()) {
                if (gcContent.get(key) > currentMax) {
                    currentMax = gcContent.get(key);
                    maxContentId = key;
                }
            }

            System.out.println(maxContentId + "\n" + gcContent.get(maxContentId));

            br.close();


        } catch (FileNotFoundException e) {
            System.out.println("File not found");
        } catch (IOException e) {
            System.out.println("Error reading file");
        }


    }

    public static void countPointMutations(String filename) {
        
        try {
            FileReader reader = new FileReader(filename);
            BufferedReader br = new BufferedReader(reader);

            String s = br.readLine();
            String t = br.readLine();

            int count = 0;
            for (int i = 0; i < s.length(); i++) {
                if (s.charAt(i) != t.charAt(i)) {
                    ++count;
                }
            }
            br.close();
            System.out.println(count);
            
        } catch (FileNotFoundException e) {
            System.out.println("File not found");
        } catch (IOException e) {
            System.out.println("IO Exception");
        }
    }
    
    public static void findMotif(String filename) {
        try {

            FileReader reader = new FileReader(filename);
            BufferedReader br = new BufferedReader(reader);

            String s = br.readLine();
            String t = br.readLine();

            ArrayList<Integer> startLocations = new ArrayList<>();

            for (int i = 0; i < s.length(); i++) {
                if (t.charAt(0) == s.charAt(i)) {
                   boolean validStartLocation = true;
                   int counterSubstring = 0;
                    if (!(i + t.length() > s.length())) {
                        for (int j = i; j < i + t.length(); j++) {
                                if (!(s.charAt(j) == t.charAt(counterSubstring))) {
                                    validStartLocation = false;
                                }
                                counterSubstring++;
                        } 
                        if (validStartLocation) startLocations.add(i);
                    }
                }
            }
            br.close();
            startLocations.stream().forEach(value -> System.out.print((value + 1) + " ")); 
            
        } catch (FileNotFoundException e) {
            System.out.println("File not found");
        } catch (IOException e) {
            System.out.println("IO Exception");
        }
    }

    public static void countPermutations(String filename) {

        try {

            FileReader reader = new FileReader(filename);
            BufferedReader br = new BufferedReader(reader);

            int number = Integer.parseInt(br.readLine());

            int[] numbers = new int[number];
            ArrayList<String> allPermutations = new ArrayList<>();

            for (int i = 0; i < number; i++) {
                numbers[i] = i + 1;
            }

            String currentCombination = "";
            
            allPermutations.add(currentCombination);
            currentCombination = "";

            for (int j = 0; j < numbers.length; j++) {
                currentCombination += j + 1 + " ";
                for (int h = 0; h < numbers.length; h++) {
                    if (!(j == h)) {
                        currentCombination += h + 1 + " ";
                    }
                }
                allPermutations.add(currentCombination);
                currentCombination = "";
            } 

            for (int j = numbers.length; j > 0; j--) {
                currentCombination += j + " ";
                for (int h = numbers.length; h > 0; h--) {
                    if (!(j == h)) {
                        currentCombination += h + " ";
                    }
                }
                allPermutations.add(currentCombination);
                currentCombination = "";
            } 

            System.out.print(number * 2);
            allPermutations.stream().forEach(value -> System.out.println(value));

        
            br.close();
        } catch (IOException e) {
            System.out.println("IO Exception");
        }

    }

    public static void rnaToProtein(String filename) {
        // define codon to protein mapping
        final Map<String, String> CODON_TABLE = Map.ofEntries(
            entry("UUU", "F"), entry("CUU", "L"), entry("AUU", "I"), entry("GUU", "V"),
            entry("UUC", "F"), entry("CUC", "L"), entry("AUC", "I"), entry("GUC", "V"),
            entry("UUA", "L"), entry("CUA", "L"), entry("AUA", "I"), entry("GUA", "V"),
            entry("UUG", "L"), entry("CUG", "L"), entry("AUG", "M"), entry("GUG", "V"),
            entry("UCU", "S"), entry("CCU", "P"), entry("ACU", "T"), entry("GCU", "A"),
            entry("UCC", "S"), entry("CCC", "P"), entry("ACC", "T"), entry("GCC", "A"),
            entry("UCA", "S"), entry("CCA", "P"), entry("ACA", "T"), entry("GCA", "A"),
            entry("UCG", "S"), entry("CCG", "P"), entry("ACG", "T"), entry("GCG", "A"),
            entry("UAU", "Y"), entry("CAU", "H"), entry("AAU", "N"), entry("GAU", "D"),
            entry("UAC", "Y"), entry("CAC", "H"), entry("AAC", "N"), entry("GAC", "D"),
            entry("UAA", ""), entry("CAA", "Q"), entry("AAA", "K"), entry("GAA", "E"),
            entry("UAG", ""), entry("CAG", "Q"), entry("AAG", "K"), entry("GAG", "E"),
            entry("UGU", "C"), entry("CGU", "R"), entry("AGU", "S"), entry("GGU", "G"),
            entry("UGC", "C"), entry("CGC", "R"), entry("AGC", "S"), entry("GGC", "G"),
            entry("UGA", ""), entry("CGA", "R"), entry("AGA", "R"), entry("GGA", "G"),
            entry("UGG", "W"), entry("CGG", "R"), entry("AGG", "R"), entry("GGG", "G")
        );
        
        try {
            FileReader reader = new FileReader(filename);       
            BufferedReader br = new BufferedReader(reader);

            String rnasequence = br.readLine();

            char[] arr = rnasequence.toCharArray();
            String aminoAcidSequence = "";
            if (arr.length % 3 == 0) {

                int startPoint = 0;

                while (startPoint < arr.length - 2) {

                    String codon = "";
                        
                    for(int i = startPoint; i < startPoint + 3; i++) {
                        codon += arr[i]; 
                    }
                
                    String aminoAcid = CODON_TABLE.get(codon);

                    aminoAcidSequence += aminoAcid;

                    startPoint += 3;
                }
            }
            br.close();
            System.out.println(aminoAcidSequence);
        } catch (IOException e) {
            System.err.println("IO Exception");
        }
    }

    public static void polypeptideWeight(String filename) {

        Map<Character, Double> aminoAcidMasses = new HashMap<>();

        // Populate the Map
        aminoAcidMasses.put('A', 71.03711);
        aminoAcidMasses.put('C', 103.00919);
        aminoAcidMasses.put('D', 115.02694);
        aminoAcidMasses.put('E', 129.04259);
        aminoAcidMasses.put('F', 147.06841);
        aminoAcidMasses.put('G', 57.02146);
        aminoAcidMasses.put('H', 137.05891);
        aminoAcidMasses.put('I', 113.08406);
        aminoAcidMasses.put('K', 128.09496);
        aminoAcidMasses.put('L', 113.08406);
        aminoAcidMasses.put('M', 131.04049);
        aminoAcidMasses.put('N', 114.04293);
        aminoAcidMasses.put('P', 97.05276);
        aminoAcidMasses.put('Q', 128.05858);
        aminoAcidMasses.put('R', 156.10111);
        aminoAcidMasses.put('S', 87.03203);
        aminoAcidMasses.put('T', 101.04768);
        aminoAcidMasses.put('V', 99.06841);
        aminoAcidMasses.put('W', 186.07931);
        aminoAcidMasses.put('Y', 163.06333);
        
        try {
            FileReader reader = new FileReader(filename);
            BufferedReader br = new BufferedReader(reader);
            double sumWeight = 0;
            String sequence = br.readLine();
            for (int i = 0; i < sequence.length(); i++) {
                sumWeight += aminoAcidMasses.get(sequence.charAt(i));
            }
            br.close();
            System.out.println(sumWeight);

        } catch (IOException e) {
            System.err.println("IO Exception");
        }
    }

    public static void rnaSplicing(String filename) {

        try {
            FileReader reader = new FileReader(filename);
            BufferedReader br = new BufferedReader(reader);

            String totalRNAString = "";
            ArrayList<String> intronSequences = new ArrayList<>();

            br.readLine();
            // extract the dna sequence
            String totalDNAString = "";
            String dnaLine = "";
            while ((dnaLine = br.readLine()) != null && !(dnaLine.startsWith(">"))) {
                totalDNAString += dnaLine;
            }

            totalRNAString = "";

            // convert the main dna template strand into an rna sequence
            for (int i = 0; i < totalDNAString.length(); i++) {
                if (totalDNAString.charAt(i) == 'T') {
                    totalRNAString = totalRNAString + "U";
                } else {
                    totalRNAString = totalRNAString + totalDNAString.charAt(i);
                }
            }

            String line = "";
            // reads introns from file
            while ((line = br.readLine()) != null) {
                if (!(line.startsWith(">"))) {
                    intronSequences.add(line);
                }
            }

            // converts the dna introns into rna sequences
            for (int i = 0; i < intronSequences.size(); i++) {
                String newIntron = "";
                for (int j = 0; j < intronSequences.get(i).length(); j++) {
                    if (intronSequences.get(i).charAt(j) == 'T') {
                        newIntron += "U";
                    } else {
                        newIntron += intronSequences.get(i).charAt(j);
                    }
                }
                intronSequences.set(i, newIntron);
            }
            
            // for each intron in the intron list, look for intron in total sequence
            String concatenatedExons = totalRNAString;

            // cuts the introns from the full rna sequence, leaving only the exons
            for (String intron : intronSequences) {
                //if (!(cutIntron(concatenatedExons, intron).equals("Intron not found"))) {
                    concatenatedExons = concatenatedExons.replace(intron, "");
                //}
            }

            char[] arr = concatenatedExons.toCharArray();
            StringBuilder splicedProteinSequence = new StringBuilder();

            int startPoint = 0;

            // retrieves the corresponding amino acid for each codon
            while (startPoint < arr.length - 2) {
                StringBuilder codon = new StringBuilder();

                for (int i = startPoint; i < startPoint + 3; i++) {
                    codon.append(arr[i]);
                }

                String aminoAcid = CODON_TABLE.get(codon.toString());

                if (aminoAcid != null) {
                    splicedProteinSequence.append(aminoAcid);
                } else {
                    System.out.println("Unknown codon: " + codon);
                }

                startPoint += 3;
            }
            br.close();
            System.out.println(splicedProteinSequence.toString());


        } catch (IOException e) {
            System.err.println("IO Exception");
        }


    }

}
