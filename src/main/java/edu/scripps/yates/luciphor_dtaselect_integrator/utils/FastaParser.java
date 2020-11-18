package edu.scripps.yates.luciphor_dtaselect_integrator.utils;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang.math.NumberUtils;

import gnu.trove.list.array.TIntArrayList;

public class FastaParser {

	public static List<String> getOutside(String seq) {
		int numOpen = 0;
		final List<String> ret = new ArrayList<String>();
		StringBuffer outside = new StringBuffer();
		for (int i = 0; i < seq.length(); i++) {
			final char charAt = seq.charAt(i);
			if (charAt == '[' || charAt == '(') {
				if (!"".equals(outside.toString())) {
					ret.add(outside.toString());
				}
				numOpen++;
			} else if (charAt == ')' || charAt == ']') {
				numOpen--;
				outside = new StringBuffer();
			} else {
				if (numOpen == 0) {
					outside.append(charAt);
				}
			}
		}
		if (!"".equals(outside.toString())) {
			ret.add(outside.toString());
		}
		return ret;
	}

	/**
	 * // R.LLLQQVSLPELPGEYSMK.V --> LLLQQVSLPELPGEYSMK
	 *
	 * @param seq
	 * @return
	 */
	public static String removeBeforeAfterAAs(String seq) {
		final String point = ".";
		if (seq.contains(point)) {
			final int firstPoint = seq.indexOf(point);
			final int lastPoint = seq.lastIndexOf(point);

			if (firstPoint != lastPoint) {
				// check that there are no numbers before or afer the points,
				// which would indicate a modification
				if ((firstPoint > 0 && NumberUtils.isNumber(String.valueOf(seq.charAt(firstPoint - 1))))
						|| (seq.length() < firstPoint + 1
								&& !NumberUtils.isNumber(String.valueOf(seq.charAt(firstPoint + 1))))) {
					return seq;
				}
				if ((lastPoint > 0 && NumberUtils.isNumber(String.valueOf(seq.charAt(lastPoint - 1))))
						|| (seq.length() < lastPoint + 1
								&& !NumberUtils.isNumber(String.valueOf(seq.charAt(lastPoint + 1))))) {
					return seq;
				}
				return seq.substring(firstPoint + 1, lastPoint);
			} else if (seq.startsWith(".")) {
				return seq.substring(1);
			} else if (seq.endsWith(".")) {
				return seq.substring(0, seq.length() - 1);
			}
		}
		return seq;
	}

	/**
	 * // R.LLLQQVSLPELPGEYSMK.V --> R
	 *
	 * @param seq
	 * @return
	 */
	public static String getBeforeSeq(String seq) {
		final String point = ".";
		if (seq.contains(point)) {
			final Integer firstPoint = getBeforeSeqPointIndex(seq);
			final Integer lastPoint = getAfterSeqPointIndex(seq);

			if (firstPoint != null && lastPoint != null && firstPoint != lastPoint) {
				final String substring = seq.substring(0, firstPoint);
				return substring;
			}
		}
		return null;
	}

	/**
	 * // R.LLLQQVSL(+80.123)PELPGEYSMK.V --> LLLQQVSL(+80.32)PELPGEYSMK
	 *
	 * @param seq
	 * @return
	 */
	public static String getSequenceInBetween(String seq) {
		final String point = ".";
		if (seq.contains(point)) {
			final int firstPoint = seq.indexOf(point);
			final int lastPoint = seq.lastIndexOf(point);

			if (firstPoint != lastPoint) {
				boolean cutSequence = true;
				// only if the previous and following character of the point are
				// not numbers, otherwise, the point is from a PTM mass
				int previous = firstPoint - 1;
				int following = firstPoint + 1;
				if (previous >= 0) {
					if (NumberUtils.isDigits(String.valueOf(seq.charAt(previous)))) {
						cutSequence = false;
					}
				}
				if (following < seq.length()) {
					if (NumberUtils.isDigits(String.valueOf(seq.charAt(following)))) {
						cutSequence = false;
					}
				}
				previous = lastPoint - 1;
				following = lastPoint + 1;
				if (previous >= 0) {
					if (NumberUtils.isDigits(String.valueOf(seq.charAt(previous)))) {
						cutSequence = false;
					}
				}
				if (following < seq.length()) {
					if (NumberUtils.isDigits(String.valueOf(seq.charAt(following)))) {
						cutSequence = false;
					}
				}
				if (cutSequence) {
					final String substring = seq.substring(firstPoint + 1, lastPoint);
					return substring;
				} else {
					return seq;
				}
			}
		}
		return seq;
	}

	/**
	 * // R.LLLQQVSLPELPGEYSMK.V --> V
	 *
	 * @param seq
	 * @return
	 */
	public static String getAfterSeq(String seq) {
		final String point = ".";
		if (seq.contains(point)) {
			final Integer firstPoint = getBeforeSeqPointIndex(seq);
			final Integer lastPoint = getAfterSeqPointIndex(seq);

			if (firstPoint != null && lastPoint != null && firstPoint != lastPoint) {
				final String substring = seq.substring(lastPoint + 1, seq.length());
				return substring;
			}
		}
		return null;
	}

	/**
	 * R.LLLQQVSLPELPGEYSMK.V --> 1 <br>
	 * LLLQQVSLPELPGEYSMK --> null
	 *
	 * @param sequence
	 * @return
	 */
	private static Integer getBeforeSeqPointIndex(String sequence) {
		final TIntArrayList allPositionsOf = StringUtils.allPositionsOf(sequence, ".");
		for (int i = 0; i < allPositionsOf.size(); i++) {
			final Integer index = allPositionsOf.get(i) - 1;
			if (index < sequence.length() - 1) {
				// is followed by a number?
				final char charAt = sequence.charAt(index + 1);
				if (!Character.isDigit(charAt)) {
					return index;
				}

			}
		}
		return null;
	}

	/**
	 * R.LLLQQVSLPELPGEYSMK.V --> 1 <br>
	 * LLLQQVSLPELPGEYSMK --> null
	 *
	 * @param sequence
	 * @return
	 */
	private static Integer getAfterSeqPointIndex(String sequence) {
		final TIntArrayList allPositionsOf = StringUtils.allPositionsOf(sequence, ".");
		for (int i = allPositionsOf.size() - 1; i >= 0; i--) {
			final Integer index = allPositionsOf.get(i) - 1;
			if (index < sequence.length() - 1) {
				// is followed by a number? then is a modification like (80.009)
				final char charAt = sequence.charAt(index + 1);
				if (!Character.isDigit(charAt)) {
					return index;
				}

			}
		}
		return null;
	}

	public static List<PTMInPeptide> getPTMsInPeptide(String fullSequence) {

		final ArrayList<edu.scripps.yates.luciphor_dtaselect_integrator.utils.PTMInPeptide> ptmsInPeptide = new ArrayList<PTMInPeptide>();
		final String sequence = cleanSequence(fullSequence);
		final List<StringPosition> tmp = FastaParser.getInside(fullSequence);
		for (final StringPosition stringPosition : tmp) {
			final int position = stringPosition.position;
			char aa = 0;
			// it can be 0 when it is n-terminal or length +1 when it is on c-term
			if (position > 0 && position <= sequence.length()) {
				aa = sequence.charAt(position - 1);
			} else {
//					log.debug("Modification at position " + position + ": " + getFullSequence());
			}
			Double deltaMass = null;
			try {
				deltaMass = Double.valueOf(stringPosition.string);
			} catch (final NumberFormatException e) {

			}
			final PTMInPeptide ptm = new PTMInPeptide(position, aa, sequence, deltaMass);
			ptmsInPeptide.add(ptm);
		}

		return ptmsInPeptide;
	}

	/**
	 * Gets a list of {@link StringPosition} objects with the information inside of
	 * parenthesis or braquets. The information is the text and the position in the
	 * text, not counting the text inside the parentheis or braquets by itself.<br>
	 * The position is based on 1, that is, starting from 1 in the first character.
	 *
	 * @param seq
	 * @return
	 */
	public static List<StringPosition> getInside(String seq) {
		int numOpen = 0;
		final List<StringPosition> ret = new ArrayList<StringPosition>();
		StringBuffer inside = new StringBuffer();
		int lastNormal = 0;
		int lenthInsides = 0;
		for (int i = 0; i < seq.length(); i++) {
			final char charAt = seq.charAt(i);
			if (charAt == '[' || charAt == '(') {
				lastNormal = i - 1 - lenthInsides;
				inside = new StringBuffer();
				numOpen++;
				lenthInsides++;
			} else if (charAt == ']' || charAt == ')') {
				numOpen--;
				if (!"".equals(inside.toString())) {
					ret.add(new StringPosition(inside.toString(), lastNormal + 1));
				}
				lenthInsides++;
			} else {
				if (numOpen > 0) {
					inside.append(charAt);
					lenthInsides++;
				}
			}
		}
		return ret;
	}

	/**
	 * This function allow to get the peptide sequence as <br>
	 * <ul>
	 * <li>K.VDLSFSPSQSLPASHAHLR.V -> VDLSFSPSQSLPASHAHLR</li>
	 * <li>R.LLLQQVSLPELPGEYSMK.V + Oxidation (M) -> LLLQQVSLPELPGEYSMK</li>
	 * <li>(-)TVAAPSVFIFPPSDEQLK(S) -> TVAAPSVFIFPPSDEQLK</li>
	 * <li>K.EKS[167.00]KESAIASTEVK.L -> EKSKESAIASTEVK</li>
	 * </ul>
	 * getting just the sequence without modifications and between the pre and post
	 * AA if available
	 *
	 * @param seq
	 * @return
	 */
	public static String cleanSequence(String seq) {
		if (seq == null)
			return null;
		final String seqTmp = seq.trim();

		if (somethingExtrangeInSequence(seqTmp)) {

			// parenthesis or brackets
			final List<String> outside = getOutside(seqTmp);
			if (!outside.isEmpty()) {
				final String tmp = appendList(outside);
				final String removeBeforeAfterAAs = removeBeforeAfterAAs(tmp);
				if (!removeBeforeAfterAAs.equals(seq)) {
					return cleanSequence(removeBeforeAfterAAs);
				}
			}

		}
		final String errorMessage = "Peptide sequence '" + seq
				+ "' is not supported. Either having not recognizable characteres or in lower case? Has it a non standard PTM enconded on it? PTMs can be encoded as in PEPTID[+45.92]E";
		AssignMass.getInstance(true);
		for (int index = 0; index < seq.length(); index++) {
			final char aa = seq.charAt(index);
			if (!AssignMass.containsMass(aa)) {
				throw new IllegalArgumentException("'" + aa + "' not recognized. " + errorMessage);
			}
		}
		if (!seq.toUpperCase().equals(seq)) {
			// it has something in lower case

			throw new IllegalArgumentException(errorMessage);

		}
		return seq.toUpperCase();

	}

	private static String appendList(List<String> list) {
		final StringBuffer sb = new StringBuffer();
		for (final String string : list) {
			sb.append(string);
		}
		return sb.toString();
	}

	public static boolean somethingExtrangeInSequence(String seq) {
		for (int i = 0; i < seq.length(); i++) {
			final char charAt = seq.charAt(i);
			if (charAt == '[' || charAt == ']')
				return true;
			if (charAt == '(' || charAt == ')')
				return true;
			if (charAt == '.')
				return true;
		}
		return false;
	}

}
