package org.metagene.genestrip.smoothy;


import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.finertree.FinerTreeMain;

import java.io.File;
import java.util.Properties;

/**
 * Command line entry point for the smoothy extension. Behaves like {@link FinerTreeMain} but builds a
 * {@link SmoothyMaker} so the finer-tree and smoothy-specific goals are available.
 *
 * @param <P> the concrete smoothy project type
 */
public abstract class SmoothyMain<P extends SmoothyProject> extends FinerTreeMain<P> {
    /**
     * Creates the smoothy command-line launcher.
     */
    protected SmoothyMain() {
    }

    /**
     * Creates the {@link SmoothyMaker} that provides the finer-tree and smoothy-specific goals for the given project.
     *
     * @param project the smoothy project
     * @return the smoothy maker for {@code project}
     */
    @Override
    protected SmoothyMaker<P> createMaker(P project) {
        return new SmoothyMaker<>(project);
    }


    /**
     * Runs the smoothy extension on the default {@link SmoothyProject} type from the command line.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        new SmoothyMain<>() {
            @Override
            protected SmoothyProject createProject(GSCommon config, String name, String key, String[] fastqFiles, String csvFile, File csvDir, File fastqResDir, String taxIds, Properties commandLineProps, GSGoalKey forGoal, String dbPath, boolean quietInit) {
                return new SmoothyProject(config, name, key, fastqFiles, csvFile, csvDir, fastqResDir, taxIds, commandLineProps, forGoal, dbPath, quietInit);
            }
        }.parseAndRun(args);
    }
}