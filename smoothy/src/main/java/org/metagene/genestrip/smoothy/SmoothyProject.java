package org.metagene.genestrip.smoothy;

import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;

import java.io.File;
import java.util.Properties;

public class SmoothyProject extends GSProject {
    /**
     * Creates a fully configured smoothy project, delegating to the {@link GSProject} constructor.
     *
     * @param config           the shared system configuration
     * @param name             the project name (also the project directory name)
     * @param key              key used as a prefix for result file names
     * @param fastqFiles       fastq/fasta input paths or URLs, may be {@code null}
     * @param csvFile          mapping file listing fastq/fasta inputs, may be {@code null}
     * @param csvDir           output directory for CSV/result files; defaults to the project's csv folder
     * @param fastqResDir      output directory for filtered fastq files, may be {@code null}
     * @param taxids           comma-separated tax ids, applied as the {@code taxids} config value
     * @param commandLineProps configuration properties supplied on the command line, may be {@code null}
     * @param forGoal          the goal the project is being set up for, used when validating config keys
     * @param dbPath           explicit database path for project-less operation, may be {@code null}
     * @param quietInit        if {@code true}, suppresses informational logging during initialization
     */
    public SmoothyProject(GSCommon config, String name, String key, String[] fastqFiles, String csvFile, File csvDir, File fastqResDir, String taxids, Properties commandLineProps, GSGoalKey forGoal, String dbPath, boolean quietInit) {
        super(config, name, key, fastqFiles, csvFile, csvDir, fastqResDir, taxids, commandLineProps, forGoal, dbPath, quietInit);
    }
}
