  600 //=========================================================================
  601 //  Fit the track, return OK if fit sucecssfull
  602 //========================================================================= 
  605   for ( int loop = 0; 3 > loop ; ++loop ) {
  606     //== Fit a parabola
  607     float s0   = 0.;
  608     float sz   = 0.;
  609     float sz2  = 0.;
  610     float sz3  = 0.;
  611     float sz4  = 0.;
  612     float sd   = 0.;
  613     float sdz  = 0.;
  614     float sdz2 = 0.;
  615     float sdz3 = 0.;
  616 
  617     float t0  = 0.;
  618     float tz  = 0.;
  619     float tz2 = 0.;
  620     float td  = 0.;
  621     float tdz = 0.;
  622 
  623     for ( PrHits::iterator itH = track.hits().begin(); track.hits().end() != itH; ++itH ) {
  624       float w = (*itH)->w();
  625       float z = (*itH)->z() - m_zreference;
  626       if ( (*itH)->dxDy() != 0 ) {
  627         if ( 0 == loop ) continue;
  628         float dy = track.deltaY( *itH );
  629         t0   += w;
  630         tz   += w * z;
  631         tz2  += w * z * z;
  632         td   += w * dy;
  633         tdz  += w * dy * z;
  634       }
  635       float d = track.distance( *itH );
  636       s0   += w;
  637       sz   += w * z;
  638       sz2  += w * z * z;
  639       sz3  += w * z * z * z;
  640       sz4  += w * z * z * z * z;
  641       sd   += w * d;
  642       sdz  += w * d * z;
  643       sdz2 += w * d * z * z;
  644       sdz3 += w * d * z * z;
  645     }
  646     float b1 = sz  * sz  - s0  * sz2;
  647     float c1 = sz2 * sz  - s0  * sz3;
  648     float d1 = sd  * sz  - s0  * sdz;
  649     float b2 = sz2 * sz2 - sz * sz3;
  650     float c2 = sz3 * sz2 - sz * sz4;
  651     float d2 = sdz * sz2 - sz * sdz2;
  652 
  653     float den = (b1 * c2 - b2 * c1 );
  654     if( fabs(den) < 1e-9 ) return false;
  655     float db  = (d1 * c2 - d2 * c1 ) / den;
  656     float dc  = (d2 * b1 - d1 * b2 ) / den;
  657     float da  = ( sd - db * sz - dc * sz2 ) / s0;
  658 
  659     float day = 0.;
  660     float dby = 0.;
  661     if ( t0 > 0. ) {
  662       float deny = (tz  * tz - t0 * tz2);
  663       day = -(tdz * tz - td * tz2) / deny;
  664       dby = -(td  * tz - t0 * tdz) / deny;
  665     }
  666 
  667     track.updateParameters( da, db, dc, day, dby );
  668     float maxChi2 = 0.;
  669     for ( PrHits::iterator itH = track.hits().begin(); track.hits().end() != itH; ++itH ) {
  670       float chi2 = track.chi2( *itH );
  671       if ( chi2 > maxChi2 ) {
  672         maxChi2 = chi2;
  673       }
  674     }
  675     if ( m_maxChi2InTrack > maxChi2 ) return true;
  676   }
  677   return false;
  678 }
  679 
