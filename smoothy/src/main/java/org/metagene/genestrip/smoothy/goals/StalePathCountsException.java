package org.metagene.genestrip.smoothy.goals;

import java.io.File;

/**
 * Thrown when stored path counts belong to another database than the one they are about to be used
 * with - typically because the database was rebuilt after the counts were written. Used against it,
 * the counts would index the nodes of a tree they were not measured on.
 * <p>
 * An {@link IllegalStateException}, so that code catching that keeps working; the fields let a caller
 * or a test tell exactly what did not match without reading the message.
 */
public class StalePathCountsException extends IllegalStateException {
    private static final long serialVersionUID = 1L;

    private final File countsFile;
    private final String countsMD5;
    private final String dbMD5;

    /**
     * Creates the exception.
     *
     * @param countsFile the file holding the stale counts
     * @param dbFile     the database file they were checked against
     * @param countsMD5  the MD5 of the database the counts were measured on
     * @param dbMD5      the MD5 of the database file at hand
     */
    public StalePathCountsException(File countsFile, File dbFile, String countsMD5, String dbMD5) {
        super("The counts in " + countsFile + " were measured on another database than " + dbFile + " (MD5 "
                + countsMD5 + " instead of " + dbMD5 + "), which has been rebuilt since. Delete " + countsFile
                + " so that the counts are measured anew.");
        this.countsFile = countsFile;
        this.countsMD5 = countsMD5;
        this.dbMD5 = dbMD5;
    }

    /**
     * Returns the file holding the stale counts.
     *
     * @return the counts file
     */
    public File getCountsFile() {
        return countsFile;
    }

    /**
     * Returns the MD5 of the database the counts were measured on.
     *
     * @return the MD5 recorded with the counts
     */
    public String getCountsMD5() {
        return countsMD5;
    }

    /**
     * Returns the MD5 of the database the counts were checked against.
     *
     * @return the database's MD5
     */
    public String getDbMD5() {
        return dbMD5;
    }
}
