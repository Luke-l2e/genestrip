package org.metagene.genestrip.smoothy.goals;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.smoothy.SmoothyProject;
import org.metagene.genestrip.store.Database;

import java.io.File;
import java.io.IOException;
import java.util.Map;

/**
 * Supplies the per-genome path counts to the goals that use them, without counting again when that has
 * been done before.
 * <p>
 * The same pattern as {@code LoadDBGoal} and the finer tree's {@code LoadKMerIndexGoal}: if
 * {@link KMerGenomePathCountsGoal} has just run, its result is taken as it is; otherwise the counts are
 * read from the file {@link KMerGenomePathCountsSERGoal} wrote, which is made first if it does not exist
 * yet. Counting reads every genome of the release in scope, while reading the file takes a moment, and
 * the counts only change when the database does - so a goal working per sample depends on this one and
 * not on the counting itself.
 * <p>
 * Before the file's counts are used, the MD5 it recorded is compared with the one of the database file:
 * counts measured on another database - one rebuilt since, say - would be read against the wrong nodes,
 * so they are refused rather than used.
 *
 * @param <P> the project type
 */
public class LoadKMerGenomePathCountsGoal<P extends SmoothyProject> extends ObjectGoal<Map<String, long[]>, P> {
    private final ObjectGoal<Map<String, long[]>, P> pathCountsGoal;
    private final File serFile;
    private final File dbFile;

    /**
     * Creates the goal.
     *
     * @param project        the project this goal belongs to
     * @param key            the key identifying this goal
     * @param pathCountsGoal the goal counting the path counts, whose result is taken if it is at hand
     * @param serGoal        the goal writing the path counts to the file read otherwise
     * @param storeDBGoal    the goal storing the database the counts belong to; only its file is used, to
     *                       check the counts against it, and it is not a dependency
     * @param deps           any further goals this goal depends on
     */
    @SafeVarargs
    public LoadKMerGenomePathCountsGoal(P project, GoalKey key, ObjectGoal<Map<String, long[]>, P> pathCountsGoal,
                                        FileGoal<P> serGoal, FileGoal<P> storeDBGoal, Goal<P>... deps) {
        super(project, key, Goal.append(deps, pathCountsGoal, serGoal));
        this.pathCountsGoal = pathCountsGoal;
        this.serFile = serGoal.getFile();
        this.dbFile = storeDBGoal.getFile();
    }

    @Override
    protected void doMakeThis() {
        if (pathCountsGoal.isMade()) {
            set(pathCountsGoal.get());
            return;
        }
        try {
            KMerGenomePathCountsSERGoal.StoredPathCounts stored = KMerGenomePathCountsSERGoal.load(serFile);
            checkDatabase(stored.getDbMD5());
            set(stored.getCounts());
        } catch (IOException | ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Refuses counts that were measured on another database than the one in the project's {@code db}
     * folder. Only the database's recorded configuration is read for that, not the database itself.
     *
     * @param countsMD5 the MD5 of the database the counts were measured on, as the file recorded it
     * @throws IOException if the database's configuration cannot be read
     */
    private void checkDatabase(String countsMD5) throws IOException {
        if (!dbFile.exists()) {
            if (getLogger().isWarnEnabled()) {
                getLogger().warn("Database " + dbFile + " not found, so the counts in " + serFile
                        + " cannot be checked against it.");
            }
            return;
        }
        String dbMD5 = Database.loadConfigInfo(dbFile).getProperty(GSProject.DB_MD5);
        if (countsMD5 == null || dbMD5 == null) {
            if (getLogger().isWarnEnabled()) {
                getLogger().warn("No database MD5 recorded, so the counts in " + serFile
                        + " cannot be checked against " + dbFile + ".");
            }
            return;
        }
        if (!countsMD5.equals(dbMD5)) {
            throw new StalePathCountsException(serFile, dbFile, countsMD5, dbMD5);
        }
    }
}
