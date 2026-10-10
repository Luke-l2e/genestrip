package org.metagene.genestrip.smoothy.goals;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.smoothy.SmoothyProject;
import org.metagene.genestrip.store.Database;

import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Writes the per-genome path counts of {@link KMerGenomePathCountsGoal} to a serialized Java object
 * file, so that later goals or analyses can reuse them without another pass over the RefSeq release.
 * {@link LoadKMerGenomePathCountsGoal} reads them back.
 * <p>
 * The file holds two objects: the MD5 of the database the counts were measured on, as the database
 * records it under {@link GSProject#DB_MD5}, and then the goal's result as a
 * {@code HashMap<String, long[]>} keyed by a genome's tax id, whose array holds the genome's k-mer
 * occurrences per position on its path to the root, {@code 0} being the genome's own node. The nodes
 * themselves are not stored, since they follow from the tax tree of that database - which is why the
 * MD5 is kept: counts read against a rebuilt database would index the wrong nodes.
 *
 * @param <P> the project type
 */
public class KMerGenomePathCountsSERGoal<P extends SmoothyProject> extends FileGoal<P> {
    private final ObjectGoal<Map<String, long[]>, P> pathCountsGoal;
    private final ObjectGoal<Database, P> dbGoal;

    /**
     * Creates the goal writing the path counts SER file.
     *
     * @param project        the project this goal belongs to
     * @param key            the key identifying this goal
     * @param pathCountsGoal the goal providing the per-genome path counts
     * @param dbGoal         the goal supplying the database the counts were measured on, for its MD5
     * @param deps           any further goals this goal depends on
     */
    @SafeVarargs
    public KMerGenomePathCountsSERGoal(P project, GoalKey key, ObjectGoal<Map<String, long[]>, P> pathCountsGoal,
                                       ObjectGoal<Database, P> dbGoal, Goal<P>... deps) {
        super(project, key, Goal.append(deps, pathCountsGoal, dbGoal));
        this.pathCountsGoal = pathCountsGoal;
        this.dbGoal = dbGoal;
    }

    /**
     * The file sits in the project's {@code db} folder beside the database it was measured on, as the
     * CSV file of {@link KMerGenomePathCountsCSVGoal} does. Its name is derived from the project and
     * goal names, as for every other output file.
     */
    @Override
    public List<File> getFiles() {
        return Collections.singletonList(getProject().getOutputFile(getProject().getDBDir(), getKey().getName(),
                null, null, GSProject.GSFileType.SER, false));
    }

    @Override
    protected void makeFile(File file) throws IOException {
        // Copied into a HashMap so that what is written is serializable whatever map the goal returns.
        HashMap<String, long[]> counts = new HashMap<>(pathCountsGoal.get());
        String dbMD5 = dbGoal.get().getConfigInfo().getProperty(GSProject.DB_MD5);
        try (ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(file))) {
            oOut.writeObject(dbMD5);
            oOut.writeObject(counts);
        }
    }

    /**
     * Reads back path counts written by this goal.
     *
     * @param file the SER file to read
     * @return the counts together with the MD5 of the database they were measured on
     * @throws IOException            if the file cannot be read, or was written by an earlier version
     *                                of this goal that did not record the database's MD5
     * @throws ClassNotFoundException if the serialized map cannot be resolved
     */
    @SuppressWarnings("unchecked")
    public static StoredPathCounts load(File file) throws IOException, ClassNotFoundException {
        try (ObjectInputStream oIn = new ObjectInputStream(StreamProvider.getInputStreamForFile(file))) {
            Object first = oIn.readObject();
            if (first != null && !(first instanceof String)) {
                throw new IOException(file + " was written by an earlier version without the database's MD5."
                        + " Delete it so that it is written anew.");
            }
            return new StoredPathCounts((String) first, (Map<String, long[]>) oIn.readObject());
        }
    }

    /**
     * Path counts as read back from a SER file, with the MD5 of the database they were measured on.
     */
    public static class StoredPathCounts {
        private final String dbMD5;
        private final Map<String, long[]> counts;

        /**
         * Creates the holder.
         *
         * @param dbMD5  the database's MD5, or {@code null} if the database recorded none
         * @param counts the per-genome path counts, keyed by the genome's tax id
         */
        public StoredPathCounts(String dbMD5, Map<String, long[]> counts) {
            this.dbMD5 = dbMD5;
            this.counts = counts;
        }

        /**
         * Returns the MD5 of the database the counts were measured on.
         *
         * @return the MD5, or {@code null} if the database recorded none
         */
        public String getDbMD5() {
            return dbMD5;
        }

        /**
         * Returns the per-genome path counts.
         *
         * @return the counts, keyed by the genome's tax id
         */
        public Map<String, long[]> getCounts() {
            return counts;
        }
    }
}
