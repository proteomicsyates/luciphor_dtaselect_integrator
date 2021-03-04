package edu.scripps.yates.luciphor_dtaselect_integrator;

import java.io.File;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.springframework.boot.CommandLineRunner;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;

@SpringBootApplication
public class LuciphorDtaselectIntegratorApplication implements CommandLineRunner {

	private static Options options;
	private static final String OPTION_LUC = "luc";
	private static final String OPTION_LFLR = "lflr";
	private static final String OPTION_GFLR = "gflr";
	private static final String OPTION_DTA = "dta";
	public static final String OPTION_REMOVE = "rem";

	public static void main(String[] args) {
		options = defineCommandLineOptions();
		try {
			if (hasOption(args, "--help") || hasOption(args, "-h")) {
				final HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp("java -jar luciphor_dtaselect_integrator-0.0.1-SNAPSHOT.jar", options);
				return;
			}
			SpringApplication.run(LuciphorDtaselectIntegratorApplication.class, args);
		} catch (final IllegalStateException e) {
			if (e.getCause() != null)
				if (e.getCause() instanceof ParseException) {
					System.err.println(e.getMessage() + ": " + e.getCause().getMessage() + "\n");
					final HelpFormatter formatter = new HelpFormatter();
					formatter.printHelp("java -jar luciphor_dtaselect_integrator-0.0.1-SNAPSHOT.jar", options);
					return;
				} else {
					e.getCause().printStackTrace();
				}
		}
	}

	private static boolean hasOption(String[] args, String option) {
		for (final String string2 : args) {
			if (string2.equals(option)) {
				return true;
			}
		}
		return false;
	}

	private static Options defineCommandLineOptions() {
		final Options options = new Options();

		// luciphor results file
		options.addOption(
				Option.builder(OPTION_LUC).desc("Full path to luciphor results file").hasArg().required().build());
		// local FLR threshold
		options.addOption(Option.builder(OPTION_LFLR).desc("Local-FLR threshold (real number from 0 to 1.0)").hasArg()
				.required(false).build());
		// global FLR threshold
		options.addOption(Option.builder(OPTION_GFLR).desc("Global-FLR threshold (real number from 0 to 1.0)").hasArg()
				.required(false).build());
		// dtaselect file
		options.addOption(
				Option.builder(OPTION_DTA).desc("Full path to dtaselect results file").hasArg().required().build());
		// remove option
		options.addOption(Option.builder(OPTION_REMOVE)
				.desc("If present, it will remove the PSMs that don't pass threshold on Luciphor's scores")
				.hasArg(false).required().build());

		return options;

	}

	@Override
	public void run(String... args) throws Exception {

		String luciphorPath = null;
		String dtaselectPath = null;
		Double lflrThreshold = null;
		Double gflrThreshold = null;
		boolean removePSMsNotPassingThreshold = false;

		final CommandLineParser parser = new DefaultParser();
		final CommandLine cmd = parser.parse(options, args);
		luciphorPath = cmd.getOptionValue(OPTION_LUC).trim();
		final File luciphorFile = new File(luciphorPath);
		if (!luciphorFile.exists()) {
			throw new ParseException("Luciphor file not found at '" + luciphorPath + "'");
		}
		dtaselectPath = cmd.getOptionValue(OPTION_DTA).trim();
		final File dtaselectFile = new File(dtaselectPath);
		if (!dtaselectFile.exists()) {
			throw new ParseException("DTASelect file not found at '" + dtaselectPath + "'");
		}

		if (cmd.hasOption(OPTION_LFLR)) {
			try {
				lflrThreshold = Double.valueOf(cmd.getOptionValue(OPTION_LFLR).trim());
				if (lflrThreshold < 0.0 || lflrThreshold > 1.0 || lflrThreshold.isNaN() || lflrThreshold.isInfinite()) {
					throw new NumberFormatException();
				}
			} catch (final NumberFormatException e) {
				throw new ParseException("Local FLR argument not valid (" + lflrThreshold
						+ "). It must be a real number between 0 and 1.0");
			}
		}

		if (cmd.hasOption(OPTION_GFLR)) {
			try {
				gflrThreshold = Double.valueOf(cmd.getOptionValue(OPTION_GFLR).trim());
				if (gflrThreshold < 0.0 || gflrThreshold > 1.0 || gflrThreshold.isNaN() || gflrThreshold.isInfinite()) {
					throw new NumberFormatException();
				}
			} catch (final NumberFormatException e) {
				throw new ParseException("Global FLR argument not valid (" + gflrThreshold
						+ "). It must be a real number between 0 and 1.0");
			}
		}
		if (cmd.hasOption(OPTION_REMOVE)) {
			removePSMsNotPassingThreshold = true;
		}
		final LuciphorDtaselectIntegrator luciphorIntegrator = new LuciphorDtaselectIntegrator(luciphorPath,
				dtaselectPath, lflrThreshold, gflrThreshold, removePSMsNotPassingThreshold);
		luciphorIntegrator.run();
		System.out.println("Program finished correctly");

	}

}
