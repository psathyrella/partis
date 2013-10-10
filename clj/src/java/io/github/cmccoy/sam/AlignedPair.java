package io.github.cmccoy.sam;

import com.google.common.base.Objects;
import com.google.common.collect.ComparisonChain;


public final class AlignedPair {
  private final int referencePosition, queryPosition;
  private final MatchesReference matchesReference;

  public AlignedPair(final int referencePosition,
                     final int queryPosition,
                     final MatchesReference matches) {
    this.referencePosition = referencePosition;
    this.queryPosition = queryPosition;
    this.matchesReference = matches;
  }

  public AlignedPair(final int referencePosition,
                     final int queryPosition) {
    this(referencePosition, queryPosition, MatchesReference.UNKNOWN);
  }

  public int getReferencePosition() { return this.referencePosition; }
  public int getQueryPosition() { return this.queryPosition; }
  public MatchesReference getMatchesReference() { return this.matchesReference; }

  @Override
  public String toString() {
    return Objects.toStringHelper(this)
        .add("reference", referencePosition)
        .add("query", queryPosition)
        .add("matches", matchesReference)
        .toString();
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(referencePosition, queryPosition);
  }

  public static enum MatchesReference {
    TRUE,
    FALSE,
    UNKNOWN
  };

  public boolean isInsertion() { return referencePosition < 0; }
  public boolean isDeletion() { return queryPosition < 0; }
  public boolean isIndel() { return isInsertion() || isDeletion(); }
  public boolean isMatch() { return matchesReference == MatchesReference.TRUE; };
  public boolean isMutation() { return matchesReference == MatchesReference.FALSE; };
  public boolean isUnknown() { return matchesReference == MatchesReference.UNKNOWN; };
}
