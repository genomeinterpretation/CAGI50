function ix = filterBasedOnReviewStatus(review)
    ix_single = contains(review, 'single');
    ix_no_conflict = contains(review, 'no_conflict');
    ix_multiple = contains( review, 'multiple');
    ix_no_assertion = contains( review, 'no_assertion');
    ix_no_interpretation = contains( review, 'no_interpretation');
    ix_expert = contains( review, 'expert');
    ix_empty =contains( review, '');
    ix = (ix_single | (ix_multiple & ix_no_conflict) | ix_expert | ix_empty)& ~(ix_no_assertion | ix_no_interpretation);
end

