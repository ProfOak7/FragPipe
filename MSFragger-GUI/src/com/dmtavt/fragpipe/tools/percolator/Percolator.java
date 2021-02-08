package com.dmtavt.fragpipe.tools.percolator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class Percolator {
    public static void main(final String[] args) {
        if (args.length == 0)
            convert_msfragger_pepXML_and_percolator_tsv_to_peptide_prophet_format(
                    Paths.get("/home/ci/percolator_test/23aug2017_hela_serum_timecourse_4mz_narrow_1.pepXML"),
                    Paths.get("/home/ci/percolator_test/psms.tsv"),
                    Paths.get("/home/ci/percolator_test/test.pep.xml"));
        else
            convert_msfragger_pepXML_and_percolator_tsv_to_peptide_prophet_format(
                    Paths.get(args[0]),
                    Paths.get(args[1]),
                    Paths.get(args[2])
            );
    }

    static String get_spectrum(final String line) {
        String spectrum = null;
        for (final String e : line.split("\\s"))
            if (e.startsWith("spectrum=")) {
                spectrum = e.substring("spectrum=\"".length(), e.length() - 1);
                break;
            }
        return spectrum.substring(0, spectrum.length() - 2);
    }

    public static void convert_msfragger_pepXML_and_percolator_tsv_to_peptide_prophet_format(
            final Path pepxml, final Path tsv, final Path output) {
        final Map<String, Double> percolator_dict = new HashMap<>();
        try (final BufferedReader brtsv = Files.newBufferedReader(tsv)) {
            final String percolator_header = brtsv.readLine();
            final List<String> colnames = Arrays.asList(percolator_header.split("\t"));
            final int indexOfPSMId = colnames.indexOf("PSMId");
            final int indexOfPEP = colnames.indexOf("posterior_error_prob");
            String line;
            while ((line = brtsv.readLine()) != null) {
                final String[] split = line.split("\t");
                final String raw_psmid = split[indexOfPSMId];
                final String psmid = raw_psmid.substring(0, raw_psmid.lastIndexOf("."));
                if (percolator_dict.containsKey(psmid))
                    throw new AssertionError();
                percolator_dict.put(psmid, Double.parseDouble(split[indexOfPEP]));
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        try (final BufferedReader brpepxml = Files.newBufferedReader(pepxml);
             final BufferedWriter out = Files.newBufferedWriter(output)) {
            String line;
            while ((line = brpepxml.readLine()) != null) {
                out.write(line + "\n");
                if (line.trim().equals("</search_summary>"))
                    break;
            }
            String spectrum;
            double one_minus_PEP = Double.NaN;
            long num_psms = 0;
            long num_psms_with_PEP = 0;
            final StringBuilder sb = new StringBuilder();
            while ((line = brpepxml.readLine()) != null) {
                if (line.trim().startsWith("<spectrum_query")) {
                    sb.setLength(0);
                    spectrum = get_spectrum(line);
                    boolean has_PEP = true;
                    try {
                        one_minus_PEP = 1 - percolator_dict.get(spectrum);
                    } catch (NullPointerException e) {
                        has_PEP = false;
                        System.out.println("spectrum = " + spectrum);
                    }
                    ++num_psms;
                    if (has_PEP) ++num_psms_with_PEP;
                    sb.append(line).append('\n');
                    while ((line = brpepxml.readLine()) != null) {
                        if (line.trim().equals("</search_hit>")) {
                            if (has_PEP)
                                sb.append(
                                        String.format(
                                                "          <analysis_result analysis=\"peptideprophet\">\n" +
                                                        "            <peptideprophet_result probability=\"%f\" all_ntt_prob=\"(%f,%f,%f)\">\n" +
                                                        "            </peptideprophet_result>\n" +
                                                        "          </analysis_result>\n", one_minus_PEP, one_minus_PEP, one_minus_PEP, one_minus_PEP));
                        }
                        sb.append(line).append("\n");
                        if (line.trim().equals("</spectrum_query>")) {
                            if (has_PEP)
                                out.write(sb.toString());
                            break;
                        }
                    }
                }
            }
            System.out.println("num_psms_with_PEP = " + num_psms_with_PEP);
            System.out.println("num_psms = " + num_psms);
            System.out.println("missing psms = " + (1 - (double) num_psms_with_PEP / num_psms));
            out.write("  </msms_run_summary>\n" +
                    "</msms_pipeline_analysis>");
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
