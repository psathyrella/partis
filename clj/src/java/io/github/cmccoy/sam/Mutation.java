package io.github.cmccoy.sam;


import com.google.common.base.Objects;
import com.google.common.collect.ComparisonChain;
import static com.google.common.base.Preconditions.*;

/**
 * An immutable mutation.
 */
public class Mutation implements Comparable<Mutation> {
    public static enum MutationType {
        INSERTION("I"), DELETION("D"), MUTATION("M");

        private final String code;
        MutationType(final String code) { this.code = code; }
        public String getCode() { return this.code; }
    }

    private final MutationType type;
    private final long position;
    private final String wt;
    private final String mut;
    public Mutation(final MutationType type,
                    final long position,
                    String wt, String mut) {
        checkArgument(position >= 0, "Negative position: %d", position);
        this.type = checkNotNull(type);
        this.position = position;
        this.wt = checkNotNull(wt);
        this.mut = checkNotNull(mut);
    }

    @Override
    public String toString() {
        return Objects.toStringHelper(this)
            .add("type", type)
            .add("position", position)
            .add("wt", wt)
            .add("mut", mut)
            .toString();
    }

    @Override
    public int hashCode() {
        return Objects.hashCode(type, position, wt, mut);
    }

    @Override
    public int compareTo(Mutation other) {
        return ComparisonChain.start()
            .compare(this.position, other.position)
            .compare(this.type, other.type)
            .compare(this.wt, other.wt)
            .compare(this.mut, other.mut)
            .result();
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) return false;
        if (obj.getClass() != getClass()) return false;
        final Mutation other = (Mutation) obj;

        return Objects.equal(other.type, this.type) &&
            other.position == this.position &&
            Objects.equal(other.wt, this.wt) &&
            Objects.equal(other.mut, this.mut);
    }

    public MutationType getType() { return type; }
    public long getPosition() { return position; }
    public String getWildType() { return wt; }
    public String getMutant() { return mut; }
}
