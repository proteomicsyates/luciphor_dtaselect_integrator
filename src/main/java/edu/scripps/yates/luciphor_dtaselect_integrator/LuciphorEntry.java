package edu.scripps.yates.luciphor_dtaselect_integrator;

import java.util.List;
import java.util.Set;

import edu.scripps.yates.luciphor_dtaselect_integrator.utils.FastaParser;
import edu.scripps.yates.luciphor_dtaselect_integrator.utils.PTMInPeptide;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.THashSet;

public class LuciphorEntry {
	private final String psmID;
	private final String predictedSequence;
	private final double localFLR;
	private final double globalFLR;
	private final double score;
	private static Set<Character> modifiedAAs = new THashSet<Character>();

	public LuciphorEntry(String psmID, String predictedSequence, double localFLR, double globalFLR, double score) {
		this.psmID = psmID;
		this.predictedSequence = predictedSequence;
		this.localFLR = localFLR;
		this.globalFLR = globalFLR;
		this.score = score;
		findModifiedAAs(predictedSequence);
	}

	private void findModifiedAAs(String predictedSequence2) {
		for (int i = 0; i < predictedSequence2.length(); i++) {
			final char aa = predictedSequence2.charAt(i);
			if (Character.isLowerCase(aa)) {
				final char upperCase = Character.toUpperCase(aa);
				if (!modifiedAAs.contains(upperCase)) {
					modifiedAAs.add(upperCase);
					System.out.println(upperCase + " is considered as modified by luciphor");
				}
			}
		}
	}

	public String getPsmID() {
		return psmID;
	}

	public String getPredictedSequence() {
		return predictedSequence;
	}

	public double getLocalFLR() {
		return localFLR;
	}

	public double getGlobalFLR() {
		return globalFLR;
	}

	/**
	 * From the predictedSequence, using the original sequence, return the predicted
	 * sequence properly formated, that is, instead of lower case aminoacids
	 * indicating that they are modified, use the actual modifications<br>
	 * From predicted sequence: QCTnVTnNITDDMRGELK, <br>
	 * and originalSequence= L.QCTN(203.079373)VTNN(203.079373)ITDDMRGELK.N, return
	 * L.QCTN(203.079373)VTN(203.079373)NITDDMRGELK.N
	 * 
	 * @param originalSequence
	 * @return
	 */
	public String getFormattedPredictedSequence(String originalSequence) {
		final StringBuilder ret = new StringBuilder();
		final String prefix = FastaParser.getBeforeSeq(originalSequence);
		final String suffix = FastaParser.getAfterSeq(originalSequence);
		originalSequence = FastaParser.getSequenceInBetween(originalSequence);
		final List<PTMInPeptide> ptms = FastaParser.getPTMsInPeptide(originalSequence);
		final TIntObjectMap<PTMInPeptide> ptmsByPosition = getPTMsByPosition(ptms);
		int indexPredictedSequence = 0;
		int numPTM = 0;
		// n-term
		if (ptmsByPosition.containsKey(0)) {
			final PTMInPeptide ptm = ptmsByPosition.get(0);
			ret.append("(" + ptm.getFormattedDeltaMass() + ")");
			numPTM++;
		}
		while (true) {
			try {
				final char aa = predictedSequence.charAt(indexPredictedSequence);
				if (ptmsByPosition.containsKey(indexPredictedSequence + 1)) {
					final PTMInPeptide ptm = ptmsByPosition.get(indexPredictedSequence + 1);
					if (!modifiedAAs.contains(Character.toUpperCase(aa))) {
						ret.append(Character.toUpperCase(aa) + "(" + ptm.getFormattedDeltaMass() + ")");
						numPTM++;
						indexPredictedSequence++;
						continue;
					}
				}

				if (Character.isLowerCase(aa)) {
					final PTMInPeptide ptm = ptms.get(numPTM++);
					ret.append(Character.toUpperCase(aa) + "(" + ptm.getFormattedDeltaMass() + ")");
				} else {
					ret.append(aa);
				}
				indexPredictedSequence++;
			} finally {
				if (indexPredictedSequence >= predictedSequence.length()) {
					break;
				}
			}
		}
		if (suffix != null) {
			ret.append(".").append(suffix);
		}
		if (prefix != null) {
			return prefix + "." + ret.toString();
		} else {
			return ret.toString();
		}
	}

	private TIntObjectMap<PTMInPeptide> getPTMsByPosition(List<PTMInPeptide> ptms) {
		final TIntObjectMap<PTMInPeptide> ret = new TIntObjectHashMap<PTMInPeptide>();
		for (final PTMInPeptide ptmInPeptide : ptms) {
			ret.put(ptmInPeptide.getPosition(), ptmInPeptide);
		}
		return ret;
	}

	public double getScore() {
		return score;
	}

}
