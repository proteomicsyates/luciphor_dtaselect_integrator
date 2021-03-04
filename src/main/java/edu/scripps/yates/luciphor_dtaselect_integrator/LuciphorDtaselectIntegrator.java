package edu.scripps.yates.luciphor_dtaselect_integrator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectIntHashMap;

public class LuciphorDtaselectIntegrator {
	private final File luciphorPath;
	private final File dtaselectPath;
	private final Double lflrThreshold;
	private final Double gflrThreshold;
	private final boolean removePSMsNotPassingThreshold;
	private boolean copied;
	private final static String COL_SPEC_ID = "specId";
	private static final String COL_PRED_PEP1 = "predictedPep1";
	private static final String COL_GLOBAL_FLR = "globalFLR";
	private static final String COL_LOCAL_FLR = "localFLR";
	private static final String DTA_COL_PSMID = "FileName";
	private static final String DTA_COL_SEQUENCE = "Sequence";
	private static final String COL_PEP1SCORE = "pep1score";
	private static final String COL_PEP2SCORE = "pep2score";
	private static final String COL_DELTASCORE = "deltaScore";
	private static final String REDUNDANCY = "Redundancy";
	private static final String SEQUENCE_COUNT = "Sequence Count";
	private static final String SPECTRUM_COUNT = "Spectrum Count";

	public LuciphorDtaselectIntegrator(String luciphorPath, String dtaselectPath, Double lflrThreshold,
			Double gflrThreshold, boolean removePSMsNotPassingThreshold) {
		this.luciphorPath = new File(luciphorPath);
		this.dtaselectPath = new File(dtaselectPath);
		this.lflrThreshold = lflrThreshold;
		this.gflrThreshold = gflrThreshold;
		this.removePSMsNotPassingThreshold = removePSMsNotPassingThreshold;
	}

	public void run() throws Exception {
		final List<LuciphorEntry> luciphorEntries = readLuciphorFile();
		final List<LuciphorEntry> filteredLuciphorEntries = filter(luciphorEntries);
		final Map<String, LuciphorEntry> filteredLuciphorEntriesByPSMID = getMapByPSMID(filteredLuciphorEntries);
		processDTASelect(filteredLuciphorEntriesByPSMID);

	}

	/**
	 * It will read the dtaselect file and will override the entries that correspond
	 * to the entries in Luciphor that pass the filters
	 * 
	 * @param filteredLuciphorEntriesByPSMID
	 * @throws IOException
	 */
	private void processDTASelect(Map<String, LuciphorEntry> filteredLuciphorEntriesByPSMID) throws IOException {
		File backupDTASelect = null;
		int numChanged = 0;
		try {
			// backup file to DTASelect-filter.txt_original
			backupDTASelect = new File(dtaselectPath.getParent() + File.separator
					+ FilenameUtils.getName(dtaselectPath.getAbsolutePath()) + "_original");
			FileUtils.copyFile(dtaselectPath, backupDTASelect);
			copied = true;
			// create the new one with the same name as the original
			final FileWriter fw = new FileWriter(this.dtaselectPath.getAbsolutePath(), false);
			// read the backup one
			final BufferedReader br = new BufferedReader(new FileReader(backupDTASelect));
			String line = null;
			boolean isPSMLine = false;
			boolean passedPSMHeader = false;
			boolean isFinalTable = false;
			boolean isProtein = false;
			TObjectIntMap<String> indexByPSMHeader = null;
			TObjectIntMap<String> indexByProteinHeader = null;
			int numPSMs = 0;
			int numPSMsRemoved = 0;
			int numProteinsRemoved = 0;
			int numPSMsOfLastProtein = 0;
			StringBuilder linesForProtein = new StringBuilder();
			while ((line = br.readLine()) != null) {

				// what is this line?
				if (isFinalTable) {
					// write it as it is
					fw.write(line + "\n");
					continue;
				} else if (line.startsWith("Locus\t")) {
					indexByProteinHeader = getIndexByHeader(line);
					// write it as it is
					fw.write(line + "\n");
					continue;
				} else if (line.startsWith("Unique\t")) {
					// get the indexes of the headers
					indexByPSMHeader = getIndexByHeader(line);
					// add new columns for the original sequence and the luciphor scores
					fw.write(line + "\toriginal_sequence\tluciphor_pep1Score\tluciphor_pep2Score\tluciphor_deltaScore\t"
							+ COL_GLOBAL_FLR + "\t" + COL_LOCAL_FLR + "\n");
					passedPSMHeader = true;
					continue;

				} else if (line.startsWith("\tProteins")) { // starts the table at the bottom of DTASelect
					if (!"".equals(linesForProtein.toString())) {
						if (numPSMsOfLastProtein > 0) {
							fw.write(processProteinLines(linesForProtein.toString(), indexByPSMHeader,
									indexByProteinHeader)); // write latest protein info
						} else {
							numProteinsRemoved++;
						}
					}
					linesForProtein = new StringBuilder();
					isPSMLine = false;
					isProtein = false;
					isFinalTable = true;
					// write it as it is
					fw.write(line + "\n");
					continue;
				} else if (line.startsWith("\t")) {
					isPSMLine = true;
					isProtein = false;
				} else if (line.startsWith("*\t")) {
					isPSMLine = true;
					isProtein = false;
				} else {
					if (passedPSMHeader && !isFinalTable) {
						// THIS IS A PROTEIN LINE

						if (!isProtein) { // if before we were not in a protein
							// this is a protein that is new
							if (!"".equals(linesForProtein.toString())) {
								if (numPSMsOfLastProtein > 0) {
									// I print the last protein lines if it had some PSMs
									fw.write(processProteinLines(linesForProtein.toString(), indexByPSMHeader,
											indexByProteinHeader));
									linesForProtein = new StringBuilder();
									numPSMsOfLastProtein = 0;
								} else {
									linesForProtein = new StringBuilder();
									numProteinsRemoved++;
								}
							}
						}
						linesForProtein.append(line + "\n");
						isProtein = true;
					} else {
						// write it as it is
						fw.write(line + "\n");
						continue;
					}
					isPSMLine = false;
				}

				if (isPSMLine) {
					numPSMs++;
					final String[] split = line.split("\t");
					final String psmID = split[indexByPSMHeader.get(DTA_COL_PSMID)];
					if (filteredLuciphorEntriesByPSMID.containsKey(psmID)) {
						final LuciphorEntry luciphorEntry = filteredLuciphorEntriesByPSMID.get(psmID);
						// replace sequence
						final String originalSequence = split[indexByPSMHeader.get(DTA_COL_SEQUENCE)];
						final String sequenceToReplace = luciphorEntry.getFormattedPredictedSequence(originalSequence);
//						if (originalSequence.equals(sequenceToReplace)) {
//							fw.write(line + "\t\t\t\t\t\t\n");
//							continue;
//						}
						// create new array of values
						final List<String> newLineValues = new ArrayList<String>();
						for (int i = 0; i < split.length; i++) {
							if (i == indexByPSMHeader.get(DTA_COL_SEQUENCE)) {
								newLineValues.add(sequenceToReplace);
							} else {
								newLineValues.add(split[i]);
							}
						}
						// now add at the end the new columns
						if (originalSequence.equals(sequenceToReplace)) {
							newLineValues.add("");
						} else {
							numChanged++;
							newLineValues.add(originalSequence);
						}

						newLineValues.add(String.valueOf(luciphorEntry.getPep1Score()));
						newLineValues.add(String.valueOf(luciphorEntry.getPep2Score()));
						newLineValues.add(String.valueOf(luciphorEntry.getDeltaScore()));
						newLineValues.add(String.valueOf(luciphorEntry.getGlobalFLR()));
						newLineValues.add(String.valueOf(luciphorEntry.getLocalFLR()));

						// get the new line as string
						final StringBuilder newLine = new StringBuilder();
						int i = 0;
						for (final String value : newLineValues) {
							if (i > 0) {
								newLine.append("\t");
							}
							newLine.append(value);
							i++;
						}

						// write the new line for the PSM
						linesForProtein.append(newLine.toString() + "\n");
						numPSMsOfLastProtein++;
						continue;
					} else {
						// in this case it is because it doesn't pass the threshold
						if (!removePSMsNotPassingThreshold) {
							linesForProtein.append(line + "\t\t\t\t\t\t\n");
							numPSMsOfLastProtein++;
						} else {
							numPSMsRemoved++;
						}
					}
				}
			}

			br.close();
			fw.close();
			System.out.println(numChanged + "/" + numPSMs + " " + getPercentageString(numChanged, numPSMs)
					+ " PSM entries with some changes in their PTM localizations were incorporated in the DTASelect file");
			if ((this.lflrThreshold != null || gflrThreshold != null) && removePSMsNotPassingThreshold) {
				System.out.println(numPSMsRemoved + " PSMs and " + numProteinsRemoved
						+ " proteins were removed because Luciphor didn't give enough confidence to them (and option "
						+ LuciphorDtaselectIntegratorApplication.OPTION_REMOVE + " was activated)");
			}
		} catch (

		final Exception e) {
			if (copied) {
				// copy back
				FileUtils.copyFile(backupDTASelect, dtaselectPath);
			}
		} finally {

		}
	}

	/**
	 * Processes the lines containing a protein or proteins of a protein group and
	 * their PSMs, modifying the columns for spec counts and peptides
	 * 
	 * @param lines
	 * @param numPSMsOfLastProtein
	 * @param sequencesPlusChargeOfLastProtein
	 * @param indexByPSMHeader
	 * @return
	 */
	private String processProteinLines(String linesString, TObjectIntMap<String> indexByPSMHeader,
			TObjectIntMap<String> indexByProteinHeader) {
		final String[] lines = linesString.split("\n");
		// first, I look for the lines that are PSMs and I grab the peptide sequences
		int numPSMEntries = 0;
		int numTotalPSMs = 0;
		for (final String line : lines) {
			if (line.startsWith("*\t") || line.startsWith("\t")) {
				// is a PSM
				numPSMEntries++;
				final String[] psmColumns = line.split("\t");
				final int spc = Integer.valueOf(psmColumns[indexByPSMHeader.get(REDUNDANCY)]);
				numTotalPSMs += spc;
			}
		}
		// now we have the number of entries and the number of PSMs
		final StringBuilder sb = new StringBuilder();
		for (final String line : lines) {
			if (line.startsWith("*\t") || line.startsWith("\t")) {
				sb.append(line + "\n");
			} else {
				// is a protein
				final String[] proteinColumns = line.split("\t");
				proteinColumns[indexByProteinHeader.get(SEQUENCE_COUNT)] = String.valueOf(numPSMEntries);
				proteinColumns[indexByProteinHeader.get(SPECTRUM_COUNT)] = String.valueOf(numTotalPSMs);
				for (int i = 0; i < proteinColumns.length; i++) {
					if (i > 0) {
						sb.append("\t");
					}
					final String proteinColumn = proteinColumns[i];
					sb.append(proteinColumn);
				}
				sb.append("\n");
			}
		}
		return sb.toString();
	}

	private static final DecimalFormat f = new DecimalFormat("#.#%");

	private String getPercentageString(int n, int total) {
		final double percentage = n * 1.0 / total;
		return "(" + f.format(percentage) + ")";
	}

	private List<LuciphorEntry> filter(List<LuciphorEntry> luciphorEntries) {
		if (this.lflrThreshold == null && gflrThreshold == null) {
			System.out.println("No threshold was defined. All PSMs from Luciphor are considered.");
			return luciphorEntries;
		}
		final int originalNumber = luciphorEntries.size();
		final Iterator<LuciphorEntry> iterator = luciphorEntries.iterator();
		while (iterator.hasNext()) {
			final LuciphorEntry luciphorentry = iterator.next();
			if (this.lflrThreshold != null) {
				if (luciphorentry.getLocalFLR() > this.lflrThreshold) {
					iterator.remove();
					continue;
				}
			}
			if (this.gflrThreshold != null) {
				if (luciphorentry.getGlobalFLR() > this.gflrThreshold) {
					iterator.remove();
					continue;
				}
			}
		}
		System.out.println(luciphorEntries.size() + "/" + originalNumber + " "
				+ getPercentageString(luciphorEntries.size(), originalNumber)
				+ " PSMs from Luciphor pass the threshold(s)");

		return luciphorEntries;
	}

	private Map<String, LuciphorEntry> getMapByPSMID(List<LuciphorEntry> list) {
		final Map<String, LuciphorEntry> ret = new THashMap<String, LuciphorEntry>();
		list.forEach(entry -> ret.put(entry.getPsmID(), entry));
		return ret;
	}

	private List<LuciphorEntry> readLuciphorFile() throws IOException {
		final List<LuciphorEntry> ret = new ArrayList<LuciphorEntry>();
		final List<String> lines = Files.readAllLines(luciphorPath.toPath());
		final String firstLine = lines.get(0);
		final TObjectIntMap<String> indexByHeader = getIndexByHeader(firstLine);
		for (int numLine = 1; numLine < lines.size(); numLine++) {
			final String line = lines.get(numLine);
			final String[] split = line.split("\t");
			if (!indexByHeader.containsKey(COL_SPEC_ID) || !indexByHeader.containsKey(COL_PRED_PEP1)) {
				continue;
			}
			final String psmID = split[indexByHeader.get(COL_SPEC_ID)];
			final String predictedSequence = split[indexByHeader.get(COL_PRED_PEP1)];
			double localFLR = Double.NaN;
			if (indexByHeader.containsKey(COL_LOCAL_FLR)) {
				localFLR = Double.valueOf(split[indexByHeader.get(COL_LOCAL_FLR)]);
			}
			double globalFLR = Double.NaN;
			if (indexByHeader.containsKey(COL_GLOBAL_FLR)) {
				globalFLR = Double.valueOf(split[indexByHeader.get(COL_GLOBAL_FLR)]);
			}
			double pep1Score = Double.NaN;
			if (indexByHeader.containsKey(COL_PEP1SCORE)) {
				pep1Score = Double.valueOf(split[indexByHeader.get(COL_PEP1SCORE)]);
			}
			double pep2Score = Double.NaN;
			if (indexByHeader.containsKey(COL_PEP2SCORE)) {
				pep2Score = Double.valueOf(split[indexByHeader.get(COL_PEP2SCORE)]);
			}
			double deltaScore = Double.NaN;
			if (indexByHeader.containsKey(COL_DELTASCORE)) {
				deltaScore = Double.valueOf(split[indexByHeader.get(COL_DELTASCORE)]);
			}
			final LuciphorEntry luciphorEntry = new LuciphorEntry(psmID, predictedSequence, localFLR, globalFLR,
					pep1Score, pep2Score, deltaScore);
			ret.add(luciphorEntry);
		}
		System.out.println(ret.size() + " PSMs read from Luciphor file");

		return ret;
	}

	private TObjectIntMap<String> getIndexByHeader(String firstLine) {
		final TObjectIntMap<String> ret = new TObjectIntHashMap<String>();
		final String[] split = firstLine.split("\t");
		for (int index = 0; index < split.length; index++) {
			ret.put(split[index], index);
		}
		return ret;
	}
}
