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
	private boolean copied;
	private final static String COL_SPEC_ID = "specId";
	private static final String COL_PRED_PEP1 = "predictedPep1";
	private static final String COL_GLOBAL_FLR = "globalFLR";
	private static final String COL_LOCAL_FLR = "localFLR";
	private static final String DTA_COL_PSMID = "FileName";
	private static final String DTA_COL_SEQUENCE = "Sequence";

	public LuciphorDtaselectIntegrator(String luciphorPath, String dtaselectPath, Double lflrThreshold,
			Double gflrThreshold) {
		this.luciphorPath = new File(luciphorPath);
		this.dtaselectPath = new File(dtaselectPath);
		this.lflrThreshold = lflrThreshold;
		this.gflrThreshold = gflrThreshold;
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
			boolean isPSMHeader = false;
			boolean passedPSMHeader = false;
			TObjectIntMap<String> indexByHeader = null;
			int numPSMs = 0;
			while ((line = br.readLine()) != null) {
				// what is this line?
				if (line.startsWith("Unique\t")) {
					isPSMHeader = true;
					passedPSMHeader = true;
				} else if (line.startsWith("\t")) {
					isPSMLine = true;
				} else if (line.startsWith("*\t")) {
					isPSMLine = true;
				} else {
					isPSMLine = false;
					isPSMHeader = false;
				}
				if (isPSMHeader) {
					// get the indexes of the headers
					indexByHeader = getIndexByHeader(line);
					// add new columns for the original sequence and the luciphor scores
					fw.write(line + "\toriginal_sequence\t" + COL_GLOBAL_FLR + "\t" + COL_LOCAL_FLR + "\n");
					continue;
				}
				if (!isPSMLine) {
					// write it as it is
					fw.write(line + "\n");
					continue;
				} else {
					numPSMs++;
					final String[] split = line.split("\t");
					final String psmID = split[indexByHeader.get(DTA_COL_PSMID)];
					if (filteredLuciphorEntriesByPSMID.containsKey(psmID)) {
						final LuciphorEntry luciphorEntry = filteredLuciphorEntriesByPSMID.get(psmID);
						// replace sequence
						final String originalSequence = split[indexByHeader.get(DTA_COL_SEQUENCE)];
						final String sequenceToReplace = luciphorEntry.getFormattedPredictedSequence(originalSequence);
//						if (originalSequence.equals(sequenceToReplace)) {
//							fw.write(line + "\t\t\t\n");
//							continue;
//						}
						// create new array of values
						final List<String> newLineValues = new ArrayList<String>();
						for (int i = 0; i < split.length; i++) {
							if (i == indexByHeader.get(DTA_COL_SEQUENCE)) {
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
						// write the new line
						fw.write(newLine.toString() + "\n");

						continue;
					} else {
						fw.write(line + "\t\t\t\n");
					}
				}
			}

			br.close();
			fw.close();
			System.out.println(numChanged + "/" + numPSMs + " " + getPercentageString(numChanged, numPSMs)
					+ " PSM entries with some changes in their PTM localizations were incorporated in the DTASelect file");
		} catch (final Exception e) {
			if (copied) {
				// copy back
				FileUtils.copyFile(backupDTASelect, dtaselectPath);
			}
		} finally {

		}
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
			final String psmID = split[indexByHeader.get(COL_SPEC_ID)];
			final String predictedSequence = split[indexByHeader.get(COL_PRED_PEP1)];
			final double localFLR = Double.valueOf(split[indexByHeader.get(COL_LOCAL_FLR)]);
			final double globalFLR = Double.valueOf(split[indexByHeader.get(COL_GLOBAL_FLR)]);
			final LuciphorEntry luciphorEntry = new LuciphorEntry(psmID, predictedSequence, localFLR, globalFLR);

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
