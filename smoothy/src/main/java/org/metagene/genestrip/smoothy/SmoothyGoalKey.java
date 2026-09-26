
package org.metagene.genestrip.smoothy;

import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.MDDescription;

import java.io.PrintStream;
import java.lang.annotation.Annotation;

/**
 * Enumeration of the goals added by the smoothy extension. Constants annotated with
 * {@link MDDescription} contribute to the generated documentation.
 */
public enum SmoothyGoalKey implements GoalKey {
    /**
     * Counts, per genome, the k-mer occurrences the database holds at every node from the genome's
     * own node up to the root.
     */
    @MDDescription("Counts, per genome, the *k*-mer occurrences the database holds at every node from the genome's own node up to the root.")
    KMER_GENOME_PATH_COUNTS("kmergenomepathcounts");

    private final boolean forUser;
    private final String name;

    SmoothyGoalKey(String name) {
        this(name, false);
    }

    SmoothyGoalKey(String name, boolean forUser) {
        this.name = name;
        this.forUser = forUser;
    }

    /**
     * Indicates whether this goal participates in transitive cleaning.
     *
     * @return {@code true}, as smoothy goals are always transitively cleaned
     */
    @Override
    public boolean isTransClean() {
        return true;
    }

    /**
     * Indicates whether this goal is meant to be invoked directly by users.
     *
     * @return whether this goal is intended to be invoked directly by users
     */
    public boolean isForUser() {
        return forUser;
    }

    /**
     * Returns the textual name of this goal as used on the command line.
     *
     * @return the goal's name
     */
    @Override
    public String getName() {
        return name;
    }

    /**
     * Returns the goal's name.
     *
     * @return the goal's name
     */
    @Override
    public String toString() {
        return name;
    }

    /**
     * Prints a Markdown table listing the smoothy goals with their user flag and descriptions.
     *
     * @param ps the stream the Markdown table is written to
     */
    public static void printGoalInfo(PrintStream ps) {
        ps.print('|');
        ps.print("Name");
        ps.print('|');
        ps.print("User Goal");
        ps.print('|');
        ps.print("Description");
        ps.print('|');
        ps.println();

        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.println();

        for (SmoothyGoalKey goalKey : SmoothyGoalKey.values()) {
            ps.print('|');
            ps.print('`');
            ps.print(goalKey.getName());
            ps.print('`');
            ps.print('|');
            ps.print(goalKey.isForUser() ? "X" : "");
            ps.print('|');
            Annotation[] annotations;
            try {
                annotations = SmoothyGoalKey.class.getField(goalKey.name()).getAnnotations();
            } catch (NoSuchFieldException | SecurityException e) {
                throw new RuntimeException(e);
            }
            for (Annotation annotation : annotations) {
                if (annotation instanceof MDDescription) {
                    ps.print(((MDDescription) annotation).value());
                    break;
                }
            }
            ps.print('|');
            ps.println();
        }
    }
}
