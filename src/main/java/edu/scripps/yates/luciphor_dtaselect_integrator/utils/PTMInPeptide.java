package edu.scripps.yates.luciphor_dtaselect_integrator.utils;

import java.text.DecimalFormat;

import org.apache.commons.lang.builder.HashCodeBuilder;

public class PTMInPeptide {

	private final Double deltaMass;
	private final static DecimalFormat formatter6Decimals = new DecimalFormat("#.######");
	private final int position;
	private final String peptideSequence;
	public final static char NULL_CHAR = '0';
	private char aa = NULL_CHAR;

	public PTMInPeptide(int position, char aa, String peptideSequence, Double deltaMass) {
		this.position = position;
		this.setAa(aa);
		this.peptideSequence = peptideSequence;
		this.deltaMass = deltaMass;
	}

	public Double getDeltaMass() {
		return deltaMass;
	}

	public String toStringExtended() {
		if (deltaMass != null) {
			final String string = getPeptideSequence() + "-" + getAa() + getPosition() + getFormattedDeltaMass();
			return string;
		} else {
			return super.toString();
		}
	}

	@Override
	public int hashCode() {
		if (deltaMass != null) {
			int hash = 23;
			hash = hash * 31 + getPosition();
			hash = hash * 31 + getPeptideSequence().hashCode();
			hash = hash * 31 + HashCodeBuilder.reflectionHashCode(deltaMass, false);
			return hash;
		}
		return super.hashCode();
	}

	public String getFormattedDeltaMass() {
		return formatter6Decimals.format(deltaMass);
	}

	public String getPeptideSequence() {
		return peptideSequence;
	}

	public int getPosition() {
		return position;
	}

	public char getAa() {
		return aa;
	}

	public void setAa(char aa) {
		this.aa = aa;
	}
}
