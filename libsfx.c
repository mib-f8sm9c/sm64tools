//A lot of this code is heavily influenced by and even directly copied from the work done by SubDrag and
// Ice Mario on the N64 Sound Tool and VADPCM decoding/encoding. This wouldn't be possible without their
// enormous contribution. Thanks guys!

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "libsfx.h"
#include "strutils.h"
#include "utils.h"

// defines

#define DEFAULT_SFX_SAMPLING_RATE 16000

// functions

static float sfx_key_table[0x100];

static wave_table * read_wave_table(unsigned char *data, unsigned int wave_offset, unsigned int sound_bank_offset)
{
   wave_table *wav = malloc(sizeof(wave_table));
   wav->unknown_1 = read_u32_be(&data[wave_offset]);
   wav->sound_offset = read_u32_be(&data[wave_offset+4]);
   
   //loop
   unsigned int loop_offset = read_u32_be(&data[wave_offset+8]);
   if(loop_offset != 0) {
     loop_offset += sound_bank_offset + 0x10;
     wav->loop = malloc(sizeof(loop_data));
     wav->loop->start = read_u32_be(&data[loop_offset]);
     wav->loop->end = read_u32_be(&data[loop_offset+4]);
     wav->loop->count = read_u32_be(&data[loop_offset+8]);
     wav->loop->unknown = read_u32_be(&data[loop_offset+12]);
     if(wav->loop->start != 0 || wav->loop->count != 0) {
       wav->loop->state = malloc(0x10 * sizeof(unsigned));
       for (int k = 0; k < 0x10; k++) {
          wav->loop->state[k] = read_u16_be(&data[loop_offset+0x10+k*2]);
       }
     }
   }
   
   //predictor
   unsigned int predictor_offset = read_u32_be(&data[wave_offset+12]);
   if(predictor_offset != 0) {
     predictor_offset += sound_bank_offset + 16;
     wav->predictor = malloc(sizeof(predictor_data));
     wav->predictor->order = read_u32_be(&data[predictor_offset]);
     wav->predictor->predictor_count = read_u32_be(&data[predictor_offset+4]);
    unsigned int num_predictor = wav->predictor->order * wav->predictor->predictor_count * 8;
     wav->predictor->data = malloc(num_predictor * sizeof(unsigned));
     for (int k = 0; k < num_predictor; k++) {
          wav->predictor->data[k] = read_u16_be(&data[predictor_offset+8+k*2]);
     }
   }
   
   wav->sound_length = read_u32_be(&data[wave_offset+16]);
   wav->unknown_2 = read_u32_be(&data[wave_offset+20]);
   
   unsigned int flt = read_u32_be(&data[wave_offset+24]);
   wav->unknown_3 = *((float*)&flt);
   wav->unknown_4 = read_u32_be(&data[wave_offset+28]);
   
   return wav;
}

void sfx_initialize_key_table()
{
   for (int x = 0; x < 0xFF; x++)
   {
      sfx_key_table[x] = (float)((60.0 - (float)x) / 12.0) * (float)((60.0 - (float)x) / 12.0); //Squared
   }
}

static unsigned char sfx_convert_ead_game_value_to_key_base(float eadKeyvalue)
{
   float keybaseReal = (((eadKeyvalue - 0.0) < 0.00001) ? 1.0f : eadKeyvalue);

   float smallestDistance = 9999999999999.0f;
   unsigned char realKey = 0;

   for (int x = 0; x < 0x100; x++)
   {
      float distance = (fabs(keybaseReal - sfx_key_table[x]));

      if (distance < smallestDistance)
      {
         smallestDistance = distance;
         realKey = x;
      }
   }

   if (realKey > 0x7F)
      realKey = 0x7F;

   return realKey;
}

// *************** //
// VADPCM Decoding //
// By Ice Mario!   //
// *************** //

static const short sfx_itable[16] =
{
   0,1,2,3,4,5,6,7,
   -8,-7,-6,-5,-4,-3,-2,-1,
};

static const short sfx_itable_half[4] =
{
	0,1,
	-2,-1,
};

static unsigned short Flip16Bit(unsigned short ShortValue)
{
	return ((ShortValue >> 8) | ((ShortValue << 8)));
}

static unsigned long Flip32Bit(unsigned long inLong)
{
	return (((inLong & 0xFF000000) >> 24) | ((inLong & 0x00FF0000) >> 8) | ((inLong & 0x0000FF00) << 8) | ((inLong & 0x000000FF) << 24));
}

unsigned long CharArrayToLong(unsigned char* currentSpot)
{
	return Flip32Bit(((unsigned long*)currentSpot)[0]);
}

static signed short sfx_sign_extend(unsigned b, // number of bits representing the number in x
                  int x      // sign extend this b-bit number to r
)
{
   

   int m = 1 << (b - 1); // mask can be pre-computed if b is fixed

   x = x & ((1 << b) - 1);  // (Skip this if bits in x above position b are already zero.)
   return (x ^ m) - m;
}

static long determineBestEncodeIndexAndPredictor_half(signed short* predictors, long numPredictors, signed short* lastSampleSet, signed short* inPCMSamples, float* bestFitIndex, long* predictorIndex)
{
        (*predictorIndex) = 0;
        float bestPredIndex = 99999999999.0f;
 
        int bestEncodeIndex = 0;
       
        for(int p = 0; p < numPredictors; p++)
        {
                signed short* tempSampleSet = malloc(8 * sizeof(signed short));
                signed short* tmp = malloc(8 * sizeof(signed short));
 
                int index = 0;
                (*bestFitIndex) = 99999999999.0f;
 
                signed short* pred1 = &predictors[p*0x10 + 0];
                signed short* pred2 = &predictors[p*0x10 + 8];
 
                for (int testIndex = 0; testIndex < 16; testIndex++)
                {
                        for (int x = 0; x < 8; x++)
                        {
                                tempSampleSet[x] = lastSampleSet[x];
                        }
 
                        float tempFitIndex = 0;
                        for (int r = 0; r < 2; r++)
                        {
                                for (int i = 0; i < 8; i++)
                                {
                                        signed short sample = inPCMSamples[(r*8)+i];
 
                                        signed long total = pred1[i] * tempSampleSet[6];
                                        total += (pred2[i] * tempSampleSet[7]);
 
                                        if (i>0)
                                        {
                                                for(int x=i-1; x>-1; x--)
                                                {
                                                        total += ( tmp[((i-1)-x)] * pred2[x] );
                                                }
                                        }
 
                                        float bestFit = 9999999999;
                                        int bestMatch = 0;
 
                                        for (int x = 0; x < 4; x++)
                                        {
                                                signed short tmpValue = (sfx_itable_half[x] << testIndex);
                                                float newValue = (((tmpValue << 0xB) + total) >> 0xB);
                                                if ((fabs((float)(sample - newValue))) < bestFit)
                                                {
                                                        bestFit = (fabs((float)(sample - newValue)));
                                                        bestMatch = x;
                                                }
                                        }
 
                                        tmp[i] = (sfx_itable_half[bestMatch] << testIndex);
                                        tempFitIndex += bestFit;
                                }
 
                       
                                for (int x = 0; x < 8; x++)
                                {
                                        tempSampleSet[x] = inPCMSamples[(r*8)+x];
                                }
                        }
 
                        if (tempFitIndex < (*bestFitIndex))
                        {
                                (*bestFitIndex) = tempFitIndex;
                                index = testIndex;
                        }
                }
 
                if ((*bestFitIndex) < bestPredIndex)
                {
                        bestPredIndex = (*bestFitIndex);
                        (*predictorIndex) = p;
                        bestEncodeIndex = index;
                }
               
                free(tmp);
                free(tempSampleSet);
        }
 
        return bestEncodeIndex;
}

static int determineBestEncodeIndexAndPredictor(signed short* predictors, int numPredictors, signed short* lastSampleSet, signed short* inPCMSamples, float* bestFitIndex, int* predictorIndex)
{
        (*predictorIndex) = 0;
        float bestPredIndex = 99999999999.0f;
 
        int bestEncodeIndex = 0;
       
        for(int p = 0; p < numPredictors; p++)
        {
                signed short* tempSampleSet = malloc(8 * sizeof(signed short));
                signed short* tmp = malloc(8 * sizeof(signed short));
 
                int index = 0;
                (*bestFitIndex) = 99999999999.0f;
 
                signed short* pred1 = &predictors[p*0x10 + 0];
                signed short* pred2 = &predictors[p*0x10 + 8];
 
                for (int testIndex = 0; testIndex < 16; testIndex++)
                {
                        for (int x = 0; x < 8; x++)
                        {
                                tempSampleSet[x] = lastSampleSet[x];
                        }
 
                        float tempFitIndex = 0;
                        for (int r = 0; r < 2; r++)
                        {
                                for (int i = 0; i < 8; i++)
                                {
                                        signed short sample = inPCMSamples[(r*8)+i];
 
                                        signed long total = pred1[i] * tempSampleSet[6];
                                        total += (pred2[i] * tempSampleSet[7]);
 
                                        if (i>0)
                                        {
                                                for(int x=i-1; x>-1; x--)
                                                {
                                                        total += ( tmp[((i-1)-x)] * pred2[x] );
                                                }
                                        }
 
                                        float bestFit = 9999999999;
                                        int bestMatch = 0;
 
                                        for (int x = 0; x < 16; x++)
                                        {
                                                signed short tmpValue = (sfx_itable[x] << testIndex);
                                                float newValue = (((tmpValue << 0xB) + total) >> 0xB);
                                                if ((fabs((float)(sample - newValue))) < bestFit)
                                                {
                                                        bestFit = (fabs((float)(sample - newValue)));
                                                        bestMatch = x;
                                                }
                                        }
 
                                        tmp[i] = (sfx_itable[bestMatch] << testIndex);
                                        tempFitIndex += bestFit;
                                }
 
                       
                                for (int x = 0; x < 8; x++)
                                {
                                        tempSampleSet[x] = inPCMSamples[(r*8)+x];
                                }
                        }
 
                        if (tempFitIndex < (*bestFitIndex))
                        {
                                (*bestFitIndex) = tempFitIndex;
                                index = testIndex;
                        }
                }
 
                if ((*bestFitIndex) < bestPredIndex)
                {
                        bestPredIndex = (*bestFitIndex);
                        (*predictorIndex) = p;
                        bestEncodeIndex = index;
                }
               
                free(tmp);
                free(tempSampleSet);
        }
 
        return bestEncodeIndex;
}

static float encode_half(signed short* inPCMSamples, int numberSamplesIn, unsigned char* outVADPCM, unsigned long* lenOut, predictor_data *predictor_info)
{
	float entropy = 0.0f;

	signed short* lastSampleSet = malloc(8 * sizeof(signed short));
	for (int x = 0; x < 8; x++)
		lastSampleSet[x] = 0x0;

	signed short* tmp = malloc(8 * sizeof(signed short));

	(*lenOut) = 0;

	for (int y = 0; y < numberSamplesIn; y += 16)
	{
		float totalBestFitDelta = 0;

		signed short* pred1;
		signed short* pred2;

		int predictor = 0;
		int index = 0;

		index = determineBestEncodeIndexAndPredictor_half(predictor_info->data, predictor_info->predictor_count, lastSampleSet, &inPCMSamples[y], &totalBestFitDelta, &predictor);

		pred1 = &predictor_info->data[predictor*0x10 + 0];
		pred2 = &predictor_info->data[predictor*0x10 + 8];

		outVADPCM[(*lenOut)++] = ((index << 4) | predictor);

		for (int r = 0; r < 2; r++)
		{
			signed short resultingValue[8];
			for (int i = 0; i < 8; i++)
			{
				signed short sample = 0;
				if ((y + (r * 8) + i) < numberSamplesIn)
                {
					sample = inPCMSamples[y+(r*8)+i];
				}

				signed long total = pred1[i] * lastSampleSet[6];
				total += (pred2[i] * lastSampleSet[7]);

				if (i>0)
				{
					for(int x=i-1; x>-1; x--)
					{
						total += ( tmp[((i-1)-x)] * pred2[x] );
					}
				}

				float bestFit = 9999999999;
				int bestMatch = 0;
				

				for (int x = 0; x < 4; x++)
				{
					signed short newValue = ((((sfx_itable_half[x] << index) << 0xB) + total) >> 0xB);
					if ((fabs((float)(sample - newValue))) < bestFit)
					{
						bestFit = (fabs((float)(sample - newValue)));
						bestMatch = x;
						resultingValue[i] = newValue;
					}
				}

				tmp[i] = (sfx_itable_half[bestMatch] << index);

				if ((i % 4) == 0)
					outVADPCM[(*lenOut)] = ((bestMatch << 6) & 0xC0);
				else if ((i % 4) == 1)
					outVADPCM[(*lenOut)] |= ((bestMatch << 4) & 0x30);
				else if ((i % 4) == 2)
					outVADPCM[(*lenOut)] |= ((bestMatch << 2) & 0x0C);
				else
				{
					outVADPCM[(*lenOut)] = (outVADPCM[(*lenOut)] | (bestMatch & 0x3));
					(*lenOut)++;
				}

				entropy += bestFit;
			}

			for (int x = 0; x < 8; x++)
			{
				//lastSampleSet[x] = inPCMSamples[y+(r*8)+x];
				lastSampleSet[x] = resultingValue[x];
			}
		}
	}


	if ((numberSamplesIn % 16) != 0)
	{
		// just cut it off for now
	}

	free(lastSampleSet);
	free(tmp);

	return entropy;
}

static float encode(signed short* inPCMSamples, int numberSamplesIn, unsigned char* outVADPCM, unsigned long* lenOut, predictor_data *predictor_info)
{
	float entropy = 0.0f;

	signed short* lastSampleSet = malloc(8 * sizeof(signed short));
	for (int x = 0; x < 8; x++)
		lastSampleSet[x] = 0x0;

	signed short* tmp = malloc(8 * sizeof(signed short));

	(*lenOut) = 0;

	for (int y = 0; y < numberSamplesIn; y += 16)
	{
		float totalBestFitDelta = 0;

		signed short* pred1;
		signed short* pred2;

		int predictor = 0;
		int index = 0;

		index = determineBestEncodeIndexAndPredictor(predictor_info->data, predictor_info->predictor_count, lastSampleSet, &inPCMSamples[y], &totalBestFitDelta, &predictor);

		pred1 = &predictor_info->data[predictor*0x10 + 0];
		pred2 = &predictor_info->data[predictor*0x10 + 8];

		//index = determineBestEncodeIndex(pred1, pred2, lastSampleSet, &inPCMSamples[y], totalBestFitDelta);

		outVADPCM[(*lenOut)++] = ((index << 4) | predictor);

		for (int r = 0; r < 2; r++)
		{
			signed short resultingValue[8];
			for (int i = 0; i < 8; i++)
			{
				signed short sample = 0;
				if ((y + (r * 8) + i) < numberSamplesIn)
                {
					sample = inPCMSamples[y+(r*8)+i];
				}

				signed long total = pred1[i] * lastSampleSet[6];
				total += (pred2[i] * lastSampleSet[7]);

				if (i>0)
				{
					for(int x=i-1; x>-1; x--)
					{
						total += ( tmp[((i-1)-x)] * pred2[x] );
					}
				}

				float bestFit = 9999999999;
				int bestMatch = 0;
				

				for (int x = 0; x < 16; x++)
				{
					signed short newValue = ((((sfx_itable[x] << index) << 0xB) + total) >> 0xB);
					if ((fabs((float)(sample - newValue))) < bestFit)
					{
						bestFit = (fabs((float)(sample - newValue)));
						bestMatch = x;
						resultingValue[i] = newValue;
					}
				}

				tmp[i] = (sfx_itable[bestMatch] << index);

				if ((i % 2) == 0)
					outVADPCM[(*lenOut)] = (bestMatch << 4);
				else
				{
					outVADPCM[(*lenOut)] = (outVADPCM[(*lenOut)] | bestMatch);
					(*lenOut)++;
				}

				entropy += bestFit;
			}

			for (int x = 0; x < 8; x++)
			{
				//lastSampleSet[x] = inPCMSamples[y+(r*8)+x];
				lastSampleSet[x] = resultingValue[x];
			}
		}
	}


	if ((numberSamplesIn % 16) != 0)
	{
		// just cut it off for now
	}

	free(lastSampleSet);
	free(tmp);

	return entropy;
}

static void decode_8( unsigned char *in, signed short *out , int index, signed short * pred1, signed short lastsmp[8])
{
   int i;
   signed short tmp[8];
   signed long total = 0;
   signed short sample =0;
   memset(out, 0, sizeof(signed short)*8);

   signed short *pred2 = (pred1 + 8);

   //printf("pred2[] = %x\n" , pred2[0]);
   for(i=0; i<8; i++)
   {
      tmp[i] = sfx_itable[(i&1) ? (*in++ & 0xf) : ((*in >> 4) & 0xf)] << index;
      tmp[i] = sfx_sign_extend(index+4, tmp[i]);
   }

   for(i=0; i<8; i++)
   {
      total = (pred1[i] * lastsmp[6]);
      total+= (pred2[i] * lastsmp[7]);

      if (i>0)
      {
         for(int x=i-1; x>-1; x--)
         {
            total += ( tmp[((i-1)-x)] * pred2[x] );
            //printf("sample: %x - pred: %x - _smp: %x\n" , ((i-1)-x) , pred2[x] , tmp[((i-1)-x)]);
         }
      }

      //printf("pred = %x | total = %x\n" , pred2[0] , total);
      float result = ((tmp[i] << 0xb) + total) >> 0xb;
      if (result > 32767)
         sample = 32767;
      else if (result < -32768)
         sample = -32768;
      else
         sample = (signed short)result;

      out[i] = sample;
   }
   // update the last sample set for subsequent iterations
   memcpy(lastsmp, out, sizeof(signed short)*8);
}

static void decode_8_half( unsigned char *in, signed short *out , int index, signed short * pred1, signed short lastsmp[8])
{
   int i;
   signed short tmp[8];
   signed long total = 0;
   signed short sample =0;
   memset(out, 0, sizeof(signed short)*8);

   signed short *pred2 = (pred1 + 8);

   //printf("pred2[] = %x\n" , pred2[0]);

   tmp[0] = (((((*in) & 0xC0) >> 6) & 0x3)) << (index);
   tmp[0] = sfx_sign_extend(index+2, tmp[0]);
   tmp[1] = (((((*in) & 0x30) >> 4) & 0x3)) << (index);
   tmp[1] = sfx_sign_extend(index+2, tmp[1]);
   tmp[2] = (((((*in) & 0x0C) >> 2) & 0x3)) << (index);
   tmp[2] = sfx_sign_extend(index+2, tmp[2]);
   tmp[3] = ((((*in++) & 0x03) & 0x3)) << (index);
   tmp[3] = sfx_sign_extend(index+2, tmp[3]);
   tmp[4] = (((((*in) & 0xC0) >> 6) & 0x3)) << (index);
   tmp[4] = sfx_sign_extend(index+2, tmp[4]);
   tmp[5] = (((((*in) & 0x30) >> 4) & 0x3)) << (index);
   tmp[5] = sfx_sign_extend(index+2, tmp[5]);
   tmp[6] = (((((*in) & 0x0C) >> 2) & 0x3)) << (index);
   tmp[6] = sfx_sign_extend(index+2, tmp[6]);
   tmp[7] = ((((*in++) & 0x03) & 0x3)) << (index);
   tmp[7] = sfx_sign_extend(index+2, tmp[7]);

   for(i=0; i<8; i++)
   {
      total = (pred1[i] * lastsmp[6]);
      total+= (pred2[i] * lastsmp[7]);

      if (i>0)
      {
         for(int x=i-1; x>-1; x--)
         {
            total += ( tmp[((i-1)-x)] * pred2[x] );
            //printf("sample: %x - pred: %x - _smp: %x\n" , ((i-1)-x) , pred2[x] , tmp[((i-1)-x)]);
         }
      }

      //printf("pred = %x | total = %x\n" , pred2[0] , total);
      float result = ((tmp[i] << 0xb) + total) >> 0xb;
      if (result > 32767)
         sample = 32767;
      else if (result < -32768)
         sample = -32768;
      else
         sample = (signed short)result;

      out[i] = sample;
   }
   // update the last sample set for subsequent iterations
   memcpy(lastsmp, out, sizeof(signed short)*8);
}

static unsigned long decode( unsigned char *in, signed short *out, unsigned long len, predictor_data *book, int decode8Only )
{
   signed short lastsmp[8];

   for (int x = 0; x < 8; x++)
      lastsmp[x] = 0;

   int index;
   int pred;
   unsigned char cmd;
   unsigned char *pin = in;
   signed short *pout = out;
   int j;
   unsigned char n1,n2;
   int total = 0;
   int _tmp;

   int samples = 0;

   // flip the predictors
   signed short *preds = (signed short*)malloc( 32 * book->predictor_count );
   for (int p = 0; p < (8 * book->order * book->predictor_count); p++)
   {
      preds[p] = book->data[p];
   }

   if (decode8Only == 0)
   {
      int _len = (len / 9) * 9;   //make sure length was actually a multiple of 9

      while(_len > 0)
      {
         index = (*in >> 4) & 0xf;
         pred = (*in & 0xf);

         // to not make zelda crash but doesn't fix it
         pred = pred % (book->predictor_count);

         _len--;

         signed short * pred1 = &preds[ pred * 16] ;

         decode_8(++in, out, index, pred1, lastsmp);
         in+=4;   _len-=4;   out+=8;

         decode_8(in, out, index, pred1, lastsmp);
         in+=4;   _len-=4;   out+=8;

         samples += 16;
      }
   }
   else
   {
      int _len = (len / 5) * 5;   //make sure length was actually a multiple of 5

      while(_len > 0)
      {
         index = (*in >> 4) & 0xf;
         pred = (*in & 0xf);

         // to not make zelda crash but doesn't fix it
         pred = pred % (book->predictor_count);

         _len--;

         signed short * pred1 = &preds[ pred * 16] ;

         decode_8_half(++in, out, index, pred1, lastsmp);
         in+=2;   _len-=2;   out+=8;

         decode_8_half(in, out, index, pred1, lastsmp);
         in+=2;   _len-=2;   out+=8;

         samples += 16;
      }
   }

   free(preds);

   return samples;
}

int extract_raw_sound_data(unsigned char* wav_file_path, unsigned char* key_base, wave_table* wav, unsigned char* final_data, unsigned long* final_count, unsigned long* sampling_rate)
{
	unsigned char* raw_data;
	int raw_length;
	unsigned long loop_start, loop_end, loop_count;
	char has_loop_data;
	
   //Open wav file path, read all data into raw_data
   FILE *file;
   
   file = fopen(wav_file_path, "rb");
   
   if (file == NULL)
   {
      //MessageBox(NULL, "Error cannot read wav file", "Error", NULL);
      fclose(file);
	  return 1;
   }
   
   fseek(file, 0, SEEK_END);
   int fileSize = ftell(file);
   rewind(file);
   
   //Check against the bad wav types
   
   unsigned char* wav_data = malloc(fileSize * sizeof(unsigned char));
   fread(wav_data, 1, fileSize, file);
   fclose(file);

   if (((((((wav_data[0] << 8) | wav_data[1]) << 8) | wav_data[2]) << 8) | wav_data[3]) != 0x52494646)
   {
      free(wav_data);
      //MessageBox(NULL, "Error not RIFF wav", "Error", NULL);
      return 2;
   }

   if (((((((wav_data[0x8] << 8) | wav_data[0x9]) << 8) | wav_data[0xA]) << 8) | wav_data[0xB]) != 0x57415645)
   {
      free(wav_data);
      //MessageBox(NULL, "Error not WAVE wav", "Error", NULL);
	  return 2;
   }

   char end_flag = 0;

   unsigned long currentOffset = 0xC;

   unsigned short channels = 0;
   (*sampling_rate) = 0;
   unsigned short bitRate = 0;

   while (end_flag == 0)
   {
      if (currentOffset >= (fileSize - 8))
         break;

      unsigned long sectionType = ((((((wav_data[currentOffset] << 8) | wav_data[currentOffset + 1]) << 8) | wav_data[currentOffset + 2]) << 8) | wav_data[currentOffset + 3]);

      if (sectionType == 0x666D7420) // fmt
      {
         unsigned long chunkSize = ((((((wav_data[currentOffset + 0x7] << 8) | wav_data[currentOffset + 0x6]) << 8) | wav_data[currentOffset + 0x5]) << 8) | wav_data[currentOffset + 0x4]);

         channels = ((wav_data[currentOffset + 0xB] << 8) | wav_data[currentOffset + 0xA]);

         if (channels != 0x0001)
         {
            //MessageBox(NULL, "Warning: Only mono wav supported", "Error", NULL);
            //end_flag = 1;
            //returnFlag = false;
            free(wav_data);
			   return 5;
         }

         (*sampling_rate) = ((((((wav_data[currentOffset + 0xF] << 8) | wav_data[currentOffset + 0xE]) << 8) | wav_data[currentOffset + 0xD]) << 8) | wav_data[currentOffset + 0xC]);
         
         bitRate = ((wav_data[currentOffset + 0x17] << 8) | wav_data[currentOffset + 0x16]);

         currentOffset += chunkSize + 8;
      }
      else if (sectionType == 0x64617461) // data
      {
         raw_length = ((((((wav_data[currentOffset + 0x7] << 8) | wav_data[currentOffset + 0x6]) << 8) | wav_data[currentOffset + 0x5]) << 8) | wav_data[currentOffset + 0x4]);

         if ((channels == 0) || ((*sampling_rate) == 0) || (bitRate == 0))
         {
            //MessageBox(NULL, "Incorrect section type (missing fmt)", "Error unknown wav format", NULL);
            free(wav_data);
			   return 3;
         }

         if (bitRate == 0x0010)
         {
            raw_data = malloc(raw_length * sizeof(unsigned char));
            for (int x = 0; x < raw_length; x++)
            {
               raw_data[x] = wav_data[currentOffset + 0x8 + x];
            }
         
            //returnFlag = true;
         }
         else
         {
            //MessageBox(NULL, "Error only 16-bit PCM wav supported", "Error", NULL);
            free(wav_data);
			   return 4;
         }

         currentOffset += raw_length + 8;
      }
      else if (sectionType == 0x736D706C) // smpl
      {
         unsigned long chunk_size = ((((((wav_data[currentOffset + 0x7] << 8) | wav_data[currentOffset + 0x6]) << 8) | wav_data[currentOffset + 0x5]) << 8) | wav_data[currentOffset + 0x4]);

         unsigned long numSampleBlocks = Flip32Bit(CharArrayToLong(&wav_data[currentOffset+0x24]));
         if (numSampleBlocks > 0)
         {
            has_loop_data = 1;

            (*key_base) = Flip32Bit(CharArrayToLong(&wav_data[currentOffset+0x14])) & 0xFF;
            loop_start = Flip32Bit(CharArrayToLong(&wav_data[currentOffset+0x34]));
            loop_end = Flip32Bit(CharArrayToLong(&wav_data[currentOffset+0x38]));
            loop_count = Flip32Bit(CharArrayToLong(&wav_data[currentOffset+0x40]));
            if (loop_count == 0)
               loop_count = 0xFFFFFFFF;
         }

         currentOffset += 8 + chunk_size;
      }
      else
      {
         unsigned long chunk_size = ((((((wav_data[currentOffset + 0x7] << 8) | wav_data[currentOffset + 0x6]) << 8) | wav_data[currentOffset + 0x5]) << 8) | wav_data[currentOffset + 0x4]);

         currentOffset += 8 + chunk_size;
      }
   }

   free(wav_data);
   
   
   ///Start C++
   
	//if (alBank->alSfx == NULL)
		//return false;

	has_loop_data = 0; //NOTE: THIS IS BLOCKING LATER CODE, TAKE OUT/RELOCATE TO START OF FUNCTION
	(*key_base) = 0x3C;
	loop_start = 0x00000000;
	loop_end = 0x00000000;
	loop_count = 0x00000000;

	wav = malloc(sizeof(wave_table));
	
	{
		//alWave->type = AL_ADPCM_WAVE;
		//alWave->adpcmWave = new ALADPCMWaveInfo();
		//alWave->adpcmWave->loop = NULL;
		/*new ALRawLoop();
		alWave->rawWave->loop->start = 0;
		alWave->rawWave->loop->end = (rawLength-2);
		alWave->rawWave->loop->count = 0;*/

		//alWave->adpcmWave->book = new ALADPCMBook();

		int numberSamples = (raw_length / 2);
		signed short* pcm_samples = malloc(numberSamples * sizeof(signed short));

		for (int x = 0; x < numberSamples; x++)
		{
			pcm_samples[x] = (signed short)(((raw_data[x*2+1] << 8)) | raw_data[x*2]);
		}

		//if (!samePred)
		{
			if (1 == 1)//decode8Only)
			{
				wav->predictor = malloc(sizeof(predictor_data));
				//alWave->adpcmWave->book->predictors = new signed short[0x10];
				wav->predictor->data = malloc(0x10 * sizeof(unsigned));
				for (int x = 0; x < 0x10; x++)
					wav->predictor->data[x] = 0x00;
					//alWave->adpcmWave->book->predictors[x] = 0x00;

				//alWave->adpcmWave->book->npredictors = 1;
				//alWave->adpcmWave->book->order = 2;
				wav->predictor->predictor_count = 1;
				wav->predictor->order = 2;
			}
			else
			{
				//TO DO: ADD THIS IN!!
				//alWave->adpcmWave->book->predictors = determineBestPredictors(alBank, alWave->adpcmWave->book->npredictors, alWave->adpcmWave->book->order, pcm_samples, numberSamples);
			}
		}
		
		//delete [] alWave->wav_data;

		unsigned long vadpcm_output_length;
		unsigned char* vadpcm_data = malloc(numberSamples * sizeof(unsigned char[numberSamples]));

		
		//TODO: ADD THESE FUNCTIONS IN
		if (1 == 1)//decode8Only)
		{
			encode_half(pcm_samples, numberSamples, vadpcm_data, vadpcm_output_length, wav->predictor);
		}
		else
		{
			encode(pcm_samples, numberSamples, vadpcm_data, vadpcm_output_length, wav->predictor);
		}

		(*final_count) = vadpcm_output_length;
		if (((*final_count) % 2) == 1)
			(*final_count)++;

		unsigned char* final_data = malloc((*final_count) * sizeof(unsigned char));
		for (int x = 0; x < (*final_count); x++)
			final_data[x] = 0x00;

		memcpy(final_data, vadpcm_data, vadpcm_output_length);

		//if ((alBank->soundBankFormat == SUPERMARIO64FORMAT)
				//|| (alBank->soundBankFormat == ZELDAFORMAT)
				//|| (alBank->soundBankFormat == STARFOX64FORMAT)
				//)
		{
			wav->loop = malloc(sizeof(loop_data));
			wav->loop->start = 0;
			wav->loop->end = ((vadpcm_output_length * 7) / 4);
			wav->loop->count = 0;

			if (has_loop_data)
			{
				wav->unknown_1 = sfx_key_table[(int)key_base];
				
			/*  Not mario format, for some reason it does this:
			
				alWave->adpcmWave->loop->count = loopCount;
				alWave->adpcmWave->loop->start = loopStart;
				alWave->adpcmWave->loop->end = loopEnd;
				alWave->adpcmWave->loop->unknown1 = 0;
				for (int x = 0; x < 0x10; x++)
					alWave->adpcmWave->loop->state[x] = alWave->adpcmWave->book->predictors[x];
			*/
			
			}
		}
		

		free(pcm_samples);
		free(vadpcm_data);
	}

	free(raw_data);

   
   //Note: there's a single return flag use above to consider: have a different return value for it??
   return 0;
}

int extract_raw_sound(unsigned char *sound_dir, unsigned char *wav_name, wave_table *wav, float key_base, unsigned char *snd_data, unsigned long sampling_rate)
{
   char bin_file[FILENAME_MAX];
   sprintf(bin_file, "%s/%s.bin", sound_dir, wav_name);
   char wav_file[FILENAME_MAX];
   sprintf(wav_file, "%s/%s.wav", sound_dir, wav_name);
   
   float sampling_rate_float = ((float)sampling_rate) / (((wav->unknown_3 - 0.0f) < 0.00001f) ? 1.0f : wav->unknown_3);

   /*if (!ignoreKeyBase)
   {
      if (
         (alBank->soundBankFormat == STANDARDFORMAT)
         || (alBank->soundBankFormat == STANDARDRNCCOMPRESSED)
         || (alBank->soundBankFormat == STANDARDRNXCOMPRESSED)
         || (alBank->soundBankFormat == BLASTCORPSZLBSTANDARD)
         || (alBank->soundBankFormat == NINDEC)
         )
      {
         samplingRateFloat = samplingRateFloat / CN64AIFCAudio::keyTable[alBank->inst[instrument]->sounds[sound]->key.keybase];
      }
      else if (
            (alBank->soundBankFormat == SUPERMARIO64FORMAT)
            || (alBank->soundBankFormat == MARIOKART64FORMAT) 
            || (alBank->soundBankFormat == ZELDAFORMAT)
            || (alBank->soundBankFormat == STARFOX64FORMAT)
         )
      {
         samplingRateFloat = samplingRateFloat / (((*reinterpret_cast<float*> (&alBank->inst[instrument]->sounds[sound]->unknown3) - 0.0) < 0.00001) ? 1.0f : *reinterpret_cast<float*> (&alBank->inst[instrument]->sounds[sound]->unknown3));
      }
   }*/

   //This algorithm is only for ADPCM WAVE format
   if ((wav == NULL) || (wav->predictor == NULL))
      return 0;
   
   unsigned char *sndData = malloc(wav->sound_length * sizeof(unsigned char));
   for(int i = 0; i < wav->sound_length; i++) {
      sndData[i] = snd_data[wav->sound_offset+i];
   }
   

   signed short* out_raw_data = malloc(wav->sound_length * 4 * sizeof(signed short));
   int n_samples = decode(sndData, out_raw_data, wav->sound_length, wav->predictor, 0);
   
   unsigned long chunk_size = 0x28 + (n_samples * 2) + 0x44 - 0x8;
   
   unsigned char *wav_data = malloc((chunk_size + 0x8 + 0x4) * sizeof(unsigned char));
   
   wav_data[0x0] = 0x52;
   wav_data[0x1] = 0x49;
   wav_data[0x2] = 0x46;
   wav_data[0x3] = 0x46;
   wav_data[0x4] = ((chunk_size >> 0) & 0xFF);
   wav_data[0x5] = ((chunk_size >> 8) & 0xFF);
   wav_data[0x6] = ((chunk_size >> 16) & 0xFF);
   wav_data[0x7] = ((chunk_size >> 24) & 0xFF);
   wav_data[0x8] = 0x57;
   wav_data[0x9] = 0x41;
   wav_data[0xA] = 0x56;
   wav_data[0xB] = 0x45;
   wav_data[0xC] = 0x66;
   wav_data[0xD] = 0x6D;
   wav_data[0xE] = 0x74;
   wav_data[0xF] = 0x20;
   wav_data[0x10] = 0x10;
   wav_data[0x11] = 0x00;
   wav_data[0x12] = 0x00;
   wav_data[0x13] = 0x00;
   wav_data[0x14] = 0x01;
   wav_data[0x15] = 0x00;
   wav_data[0x16] = 0x01;
   wav_data[0x17] = 0x00;
   wav_data[0x18] = (((unsigned long)sampling_rate_float >> 0) & 0xFF);
   wav_data[0x19] = (((unsigned long)sampling_rate_float >> 8) & 0xFF);
   wav_data[0x1A] = (((unsigned long)sampling_rate_float >> 16) & 0xFF);
   wav_data[0x1B] = (((unsigned long)sampling_rate_float >> 24) & 0xFF);
   wav_data[0x1C] = ((((unsigned long)sampling_rate_float * 2) >> 0) & 0xFF);
   wav_data[0x1D] = ((((unsigned long)sampling_rate_float * 2) >> 8) & 0xFF);
   wav_data[0x1E] = ((((unsigned long)sampling_rate_float * 2) >> 16) & 0xFF);
   wav_data[0x1F] = ((((unsigned long)sampling_rate_float * 2) >> 24) & 0xFF);
   wav_data[0x20] = 0x02;
   wav_data[0x21] = 0x00;
   wav_data[0x22] = 0x10;
   wav_data[0x23] = 0x00;
   wav_data[0x24] = 0x64;
   wav_data[0x25] = 0x61;
   wav_data[0x26] = 0x74;
   wav_data[0x27] = 0x61;

   unsigned long length = (n_samples * 2);
   
   wav_data[0x28] = ((length >> 0) & 0xFF);
   wav_data[0x29] = ((length >> 8) & 0xFF);
   wav_data[0x2A] = ((length >> 16) & 0xFF);
   wav_data[0x2B] = ((length >> 24) & 0xFF);
   
   for (int x = 0; x < n_samples; x++)
   {
      wav_data[0x2C + x * 2] = ((out_raw_data[x] >> 0) & 0xFF);
      wav_data[0x2C + x * 2 + 1] = ((out_raw_data[x] >> 8) & 0xFF);
   }

   unsigned long wavIndex = 0x2C + n_samples * 2;
   
   if (wav->loop->start != 0 || wav->loop->count != 0)
   {
      for (int x = 0; x < 0x44; x++)
         wav_data[wavIndex + x] = 0x00;

      wav_data[wavIndex + 0x0] = 0x73;
      wav_data[wavIndex + 0x1] = 0x6D;
      wav_data[wavIndex + 0x2] = 0x70;
      wav_data[wavIndex + 0x3] = 0x6C;
      wav_data[wavIndex + 0x4] = 0x3C;
      wav_data[wavIndex + 0x5] = 0x00;
      wav_data[wavIndex + 0x6] = 0x00;
      wav_data[wavIndex + 0x7] = 0x00;
      
      //This value only holds true for Mario/Zelda/StarFox formats
      {
         wav_data[wavIndex + 0x14] = sfx_convert_ead_game_value_to_key_base(key_base);
      }
      
      wav_data[wavIndex + 0x15] = 0x00;
      wav_data[wavIndex + 0x16] = 0x00;
      wav_data[wavIndex + 0x17] = 0x00;
      
      wav_data[wavIndex + 0x24] = 0x01;
      wav_data[wavIndex + 0x25] = 0x00;
      wav_data[wavIndex + 0x26] = 0x00;
      wav_data[wavIndex + 0x27] = 0x00;

      if (wav->loop->count > 0)
      {
         wav_data[wavIndex + 0x34] = ((wav->loop->start >> 0) & 0xFF);
         wav_data[wavIndex + 0x35] = ((wav->loop->start >> 8) & 0xFF);
         wav_data[wavIndex + 0x36] = ((wav->loop->start >> 16) & 0xFF);
         wav_data[wavIndex + 0x37] = ((wav->loop->start >> 24) & 0xFF);
         wav_data[wavIndex + 0x38] = (((wav->loop->end) >> 0) & 0xFF);
         wav_data[wavIndex + 0x39] = (((wav->loop->end) >> 8) & 0xFF);
         wav_data[wavIndex + 0x3A] = (((wav->loop->end) >> 16) & 0xFF);
         wav_data[wavIndex + 0x3B] = (((wav->loop->end) >> 24) & 0xFF);

         if (wav->loop->count == 0xFFFFFFFF)
         {
            wav_data[wavIndex + 0x40] = 0x00;
            wav_data[wavIndex + 0x41] = 0x00;
            wav_data[wavIndex + 0x42] = 0x00;
            wav_data[wavIndex + 0x43] = 0x00;
         }
         else
         {
            wav_data[wavIndex + 0x40] = (((wav->loop->count) >> 0) & 0xFF);
            wav_data[wavIndex + 0x41] = (((wav->loop->count) >> 8) & 0xFF);
            wav_data[wavIndex + 0x42] = (((wav->loop->count) >> 16) & 0xFF);
            wav_data[wavIndex + 0x43] = (((wav->loop->count) >> 24) & 0xFF);
         }
      }
   }
   else
   {
      for (int x = 0; x < 0x44; x++)
         wav_data[wavIndex + x] = 0x00;

      wav_data[wavIndex + 0x0] = 0x73;
      wav_data[wavIndex + 0x1] = 0x6D;
      wav_data[wavIndex + 0x2] = 0x70;
      wav_data[wavIndex + 0x3] = 0x6C;
      wav_data[wavIndex + 0x4] = 0x3C;
      wav_data[wavIndex + 0x5] = 0x00;
      wav_data[wavIndex + 0x6] = 0x00;
      wav_data[wavIndex + 0x7] = 0x00;
      
      //This value only holds true for Mario/Zelda/StarFox formats
      {
         wav_data[wavIndex + 0x14] = sfx_convert_ead_game_value_to_key_base(key_base);
      }
      wav_data[wavIndex + 0x15] = 0x00;
      wav_data[wavIndex + 0x16] = 0x00;
      wav_data[wavIndex + 0x17] = 0x00;
   }

   //Write bin file first, so file creation dates can be used to see if the wav file's been updated
   write_file(bin_file, sndData, wav->sound_length);
   
   write_file(wav_file, wav_data, (chunk_size + 0x8 + 0x4));
   
   free(wav_data);
   free(out_raw_data);

   return 1;
}

int write_sound_data(sound_data_header sound_data, unsigned char *bin_file) {
   
   unsigned long sound_data_length = 4;
   unsigned i, j;
   int ret_val = 0;
   
   sound_data_length += 8 * sound_data.data_count;
   sound_data_length += 0x10 - (sound_data_length % 0x10); //add buffer to 0x10 address
   
   for (i = 0; i < sound_data.data_count; i++) {
      sound_data_length += sound_data.data_length[i];
   }
   
   unsigned char *sndData = malloc(sound_data_length);
   memset(sndData, 0, sound_data_length);
   
   write_u16_be(&sndData[0], sound_data.unknown);
   write_u16_be(&sndData[2], sound_data.data_count);
   
   unsigned int sound_data_offset = 4 + 8 * sound_data.data_count;
   sound_data_offset += 0x10 - (sound_data_offset % 0x10); //add buffer to 0x10 address
   
   for (i = 0; i < sound_data.data_count; i++) {
      write_u32_be(&sndData[4 + i * 8], sound_data_offset);
      write_u32_be(&sndData[8 + i * 8], sound_data.data_length[i]);
      //NEED TO MAKE SURE THE DATA BLOCK ISN'T THE SAME AS A PREVIOUS ONE!!
      for(int j = 0; j < sound_data.data_length[i]; j++) {
         sndData[sound_data_offset + j] = sound_data.data[i][j];
      }
      sound_data_offset += sound_data.data_length[i];
   }
   
   //STILL NEED TO WRITE THE FILE OPEN HERE!!
   FILE *out;
   out = fopen(bin_file, "wb");
   if (out == NULL) {
      return 1;
   }
   
   int bytes_written = fwrite(sndData, 1, sound_data_length, out);
   if (bytes_written != sound_data_length) {
      ret_val = 2;
   }

   // clean up
   fclose(out);
   return ret_val;
}
 
sound_data_header read_sound_data(unsigned char *data, unsigned int data_offset) {
   
   unsigned i, j;
   sound_data_header sound_data;
   
   sound_data.unknown = read_u16_be(&data[data_offset]);
   sound_data.data_count = read_u16_be(&data[data_offset+2]);
   
   if (sound_data.data_count > 0) {
      sound_data.data = malloc(sound_data.data_count * sizeof(*sound_data.data));
      sound_data.data_length = malloc(sound_data.data_count * sizeof(unsigned int));
      for (i = 0; i < sound_data.data_count; i++) {
         unsigned int sound_data_offset = read_u32_be(&data[data_offset+i*8+4]) + data_offset;
         unsigned int sound_data_length = read_u32_be(&data[data_offset+i*8+8]);

         sound_data.data[i] = malloc(sound_data_length * sizeof(unsigned char));
         for(j = 0; j < sound_data_length; j++) {
            sound_data.data[i][j] = data[sound_data_offset+j];
         }
         sound_data.data_length[i] = sound_data_length;
      }
   }
   
   return sound_data;
}
   
int get_wav_table_size(wave_table *wav) {
   int length = 0x20;
   if(wav->loop != NULL) {
      length += 0x10;
      if(wav->loop->start != 0 || wav->loop->count != 0)
         length += 0x20;
   }
   if(wav->predictor != NULL) {
      length += 8 + wav->predictor->order * wav->predictor->predictor_count * 8 * sizeof(unsigned);
   }
   return length;
}
   
int write_sound_bank_wav(wave_table *wav, int *offset, int bankOffset, unsigned char *data) {
   int wavOffset = (*offset);
   int startWavOffset = wavOffset;
   write_u32_be(&data[wavOffset], wav->unknown_1);
   write_u32_be(&data[wavOffset + 0x4], wav->sound_offset);
   write_u32_be(&data[wavOffset + 0x10], wav->sound_length);
   write_u32_be(&data[wavOffset + 0x14], wav->unknown_2);
   write_u32_be(&data[wavOffset + 0x18], *((unsigned int*)&wav->unknown_3)); //float to uint
   write_u32_be(&data[wavOffset + 0x1C], wav->unknown_4);
   
   wavOffset += 0x20;
   
   if(wav->predictor != NULL) {
      write_u32_be(&data[startWavOffset + 0xC], wavOffset - bankOffset);
      
      write_u32_be(&data[wavOffset], wav->predictor->order);
      write_u32_be(&data[wavOffset + 0x4], wav->predictor->predictor_count);
      
      for(int i = 0; i < wav->predictor->order * wav->predictor->predictor_count * 8; i++ ) {
            write_u16_be(&data[wavOffset + 8 + i * 2], wav->predictor->data[i]);
      }
      
      wavOffset += 8 + wav->predictor->order * wav->predictor->predictor_count * 8 * 2;
      if((wavOffset % 0x10) != 0)
         wavOffset += 0x10 - (wavOffset % 0x10);
   }
   
   if(wav->loop != NULL) {
      write_u32_be(&data[startWavOffset + 0x8], wavOffset - bankOffset);
      
      write_u32_be(&data[wavOffset], wav->loop->start);
      write_u32_be(&data[wavOffset + 0x4], wav->loop->end);
      write_u32_be(&data[wavOffset + 0x8], wav->loop->count);
      write_u32_be(&data[wavOffset + 0xC], wav->loop->unknown);
      wavOffset += 0x10;
      if(wav->loop->start != 0 || wav->loop->count != 0) {
         for(int i = 0; i < 0x10; i++)
            write_u16_be(&data[wavOffset + i * 2], wav->loop->state[i]);
         wavOffset += 0x20;
      }
   }
   
   (*offset) = wavOffset;
}
   
int write_sound_bank_adrs(unsigned *adrs, int *offset, unsigned char *data) {
   write_u16_be(&data[(*offset)], adrs[0]);
   write_u16_be(&data[(*offset) + 0x2], adrs[1]);
   write_u16_be(&data[(*offset) + 0x4], adrs[2]);
   write_u16_be(&data[(*offset) + 0x6], adrs[3]);
   write_u16_be(&data[(*offset) + 0x8], adrs[4]);
   write_u16_be(&data[(*offset) + 0xA], adrs[5]);
   write_u16_be(&data[(*offset) + 0xC], adrs[6]);
   write_u16_be(&data[(*offset) + 0xE], adrs[7]);
   (*offset) += 0x10;
}
   
//TODO: FIGURE OUT A WAY TO USE DUPLICATE ADRS DATA TO SAVE ON SPACE
int write_sound_bank(sound_bank_header sound_bank, unsigned char *bin_file) {
   
   unsigned long sound_bank_length = 4;
   unsigned i, j;
   int ret_val = 0;
   
   sound_bank_length += 8 * sound_bank.bank_count;
   if((sound_bank_length % 0x10) != 0)
      sound_bank_length += 0x10 - (sound_bank_length % 0x10);
   
   for (i = 0; i < sound_bank.bank_count; i++) {
      sound_bank_length += 8; //bank list
      sound_bank_length += 14 + 4 * sound_bank.banks[i].instrument_count; //bank
      for (j = 0; j < sound_bank.banks[i].instrument_count; j++) {
         sound_bank_length += 0x20; //sound
         //need to handle identical adrs in the future here : /
         if(sound_bank.banks[i].sounds[j].adrs != NULL)
            sound_bank_length += 0x10;
         if(sound_bank.banks[i].sounds[j].wav != NULL)
            sound_bank_length += get_wav_table_size(sound_bank.banks[i].sounds[j].wav);
         if(sound_bank.banks[i].sounds[j].wav_prev != NULL)
            sound_bank_length += get_wav_table_size(sound_bank.banks[i].sounds[j].wav_prev);
         if(sound_bank.banks[i].sounds[j].wav_sec != NULL)
            sound_bank_length += get_wav_table_size(sound_bank.banks[i].sounds[j].wav_sec);
      }
      sound_bank_length += 4 * sound_bank.banks[i].percussion_count; //percussions
      for (j = 0; j < sound_bank.banks[i].percussion_count; j++) {
         sound_bank_length += 0x10;
         if(sound_bank.banks[i].percussions.items[j].adrs != NULL)
            sound_bank_length += 0x10;
         if(sound_bank.banks[i].percussions.items[j].wav != NULL)
            sound_bank_length += get_wav_table_size(sound_bank.banks[i].percussions.items[j].wav);
      }
   }
   
   //handle buffers here
   sound_bank_length += 0x1000;
   
   //Now actually write it
   unsigned char *sndData = malloc(sound_bank_length);
   memset(sndData, 0, sound_bank_length);
   
   write_u16_be(&sndData[0], sound_bank.unknown);
   write_u16_be(&sndData[2], sound_bank.bank_count);
   
   int startInstrumentOffset = 4 + sound_bank.bank_count * 8;
   if((startInstrumentOffset % 0x10) != 0)
      startInstrumentOffset += 0x10 - (startInstrumentOffset % 0x10);
   
   int instrumentOffset = 0;
   
   for (i = 0; i < sound_bank.bank_count; i++) {
      int bankOffset = startInstrumentOffset;
      int length = 0x20;
      
      //instrument info
      write_u32_be(&sndData[bankOffset], sound_bank.banks[i].instrument_count); 
      write_u32_be(&sndData[bankOffset+0x4], sound_bank.banks[i].percussion_count); 
      write_u32_be(&sndData[bankOffset+0x8], sound_bank.banks[i].unknown_1); 
      write_u32_be(&sndData[bankOffset+0xC], sound_bank.banks[i].unknown_2); 
      
      //only need to reference this one
      bankOffset += 0x10;
      
      instrumentOffset = bankOffset + (sound_bank.banks[i].instrument_count + 1) * 4;
      if((instrumentOffset % 0x10) != 0)
         instrumentOffset += 0x10 - (instrumentOffset % 0x10); //THIS NEEDS TO GO ABOVE IN THE SIZE CALCULATION
      
      for (j = 0; j < sound_bank.banks[i].instrument_count; j++) {
         //write inst data
         int headInstOffset = instrumentOffset;
         write_u32_be(&sndData[headInstOffset], sound_bank.banks[i].sounds[j].unknown);
         write_u32_be(&sndData[headInstOffset+0xC], *((unsigned int*)&sound_bank.banks[i].sounds[j].key_base_prev)); //float to uint
         write_u32_be(&sndData[headInstOffset+0x14], *((unsigned int*)&sound_bank.banks[i].sounds[j].key_base));
         write_u32_be(&sndData[headInstOffset+0x1C], *((unsigned int*)&sound_bank.banks[i].sounds[j].key_base_sec));
         
         instrumentOffset += 0x20;
         
         if(sound_bank.banks[i].sounds[j].adrs != NULL) {
            //Make this into a function please
            write_u32_be(&sndData[headInstOffset + 0x4], instrumentOffset - bankOffset);
            write_sound_bank_adrs(sound_bank.banks[i].sounds[j].adrs, &instrumentOffset, sndData);
         }
         
         if(sound_bank.banks[i].sounds[j].wav_prev != NULL) {
            write_u32_be(&sndData[headInstOffset + 0x8], instrumentOffset - bankOffset);
            write_sound_bank_wav(sound_bank.banks[i].sounds[j].wav_prev, &instrumentOffset, bankOffset, sndData);
         }
         else {
            write_u32_be(&sndData[headInstOffset + 0x8], 0);
         }
         
         if(sound_bank.banks[i].sounds[j].wav != NULL) {
            write_u32_be(&sndData[headInstOffset + 0x10], instrumentOffset - bankOffset);
            write_sound_bank_wav(sound_bank.banks[i].sounds[j].wav, &instrumentOffset, bankOffset, sndData);
         }
         
         if(sound_bank.banks[i].sounds[j].wav_sec != NULL) {
            write_u32_be(&sndData[headInstOffset + 0x18], instrumentOffset - bankOffset);
            write_sound_bank_wav(sound_bank.banks[i].sounds[j].wav_sec, &instrumentOffset, bankOffset, sndData);
         }
         
         //Write the offset in the inst table
         write_u32_be(&sndData[bankOffset + (j + 1) * 4], headInstOffset - bankOffset);
      }
      
      //percussions
      if(sound_bank.banks[i].percussion_count > 0) {
         if((instrumentOffset % 0x10) != 0)
            instrumentOffset += 0x10 - (instrumentOffset % 0x10); //THIS NEEDS TO GO ABOVE IN THE SIZE CALCULATION
         
         int percussionTableOffset = instrumentOffset;
         
         instrumentOffset += sound_bank.banks[i].percussion_count * 4;
         if((instrumentOffset % 0x10) != 0)
            instrumentOffset += 0x10 - (instrumentOffset % 0x10); //THIS NEEDS TO GO ABOVE IN THE SIZE CALCULATION
         
         for (j = 0; j < sound_bank.banks[i].percussion_count; j++) {
            //table offset
            write_u32_be(&sndData[percussionTableOffset + 4 * i], instrumentOffset - bankOffset);
            
            int headPercOffset = instrumentOffset;
            sndData[headPercOffset] = sound_bank.banks[i].percussions.items[j].unknown_1;
            sndData[headPercOffset + 1] = sound_bank.banks[i].percussions.items[j].pan;
            write_u16_be(&sndData[headPercOffset + 0x2], sound_bank.banks[i].percussions.items[j].unknown_2);
            write_u32_be(&sndData[headPercOffset + 0x8], *((unsigned int*)&sound_bank.banks[i].percussions.items[j].key_base));
            
            instrumentOffset += 0x10;
            
            if(sound_bank.banks[i].percussions.items[j].adrs != NULL) {
               write_u32_be(&sndData[headPercOffset + 0xC], instrumentOffset - bankOffset);
                write_sound_bank_adrs(sound_bank.banks[i].percussions.items[j].adrs, &instrumentOffset, sndData);
            }
            if(sound_bank.banks[i].percussions.items[j].wav != NULL) {
               write_u32_be(&sndData[headPercOffset + 0x4], instrumentOffset - bankOffset);
               write_sound_bank_wav(sound_bank.banks[i].percussions.items[j].wav, &instrumentOffset, bankOffset, sndData);
            }
            
         }
         
         write_u32_be(&sndData[bankOffset], percussionTableOffset - bankOffset);
      }
            
      if((instrumentOffset % 0x10) != 0)
         instrumentOffset += 0x10 - (instrumentOffset % 0x10); //THIS NEEDS TO GO ABOVE IN THE SIZE CALCULATION
      
      write_u32_be(&sndData[4 + i * 8], startInstrumentOffset); 
      write_u32_be(&sndData[8 + i * 8], instrumentOffset - startInstrumentOffset);
      startInstrumentOffset = instrumentOffset;
   }
   
   //STILL NEED TO WRITE THE FILE OPEN HERE!!
   FILE *out;
   out = fopen(bin_file, "wb");
   if (out == NULL) {
      return 1;
   }
   
   //by this point, instrumentOffset is a true data length value for the array
   
   int bytes_written = fwrite(sndData, 1, instrumentOffset, out);
   if (bytes_written != instrumentOffset) {
      ret_val = 2;
   }

   // clean up
   fclose(out);
   return ret_val;
}
 
sound_bank_header read_sound_bank(unsigned char *data, unsigned int data_offset) {
   
   unsigned i, j, k;
   sound_bank_header sound_banks;
   
   sound_banks.unknown = read_u16_be(&data[data_offset]);
   sound_banks.bank_count = read_u16_be(&data[data_offset+2]);
   if (sound_banks.bank_count > 0) {
      sound_banks.banks = malloc(sound_banks.bank_count * sizeof(*sound_banks.banks));
      for (i = 0; i < sound_banks.bank_count; i++) {
        unsigned int sound_bank_offset = read_u32_be(&data[data_offset+i*8+4]) + data_offset;
        //unsigned int length = read_u32_be(&data[secCtl->start+i*8+8]);
       
       sound_banks.banks[i].instrument_count = read_u32_be(&data[sound_bank_offset]);
       sound_banks.banks[i].percussion_count = read_u32_be(&data[sound_bank_offset+4]);
       sound_banks.banks[i].unknown_1 = read_u32_be(&data[sound_bank_offset+8]);
       sound_banks.banks[i].unknown_2 = read_u32_be(&data[sound_bank_offset+12]);
       
       //sounds
       if (sound_banks.banks[i].instrument_count > 0) {
          sound_banks.banks[i].sounds = malloc(sound_banks.banks[i].instrument_count * sizeof(*sound_banks.banks[i].sounds));
         for (j = 0; j < sound_banks.banks[i].instrument_count; j++) {
            unsigned int sound_offset = read_u32_be(&data[sound_bank_offset+20+j*4]);
            
            if(sound_offset != 0)
            {
               sound_offset += sound_bank_offset + 16;
              
               sound_banks.banks[i].sounds[j].unknown = read_u32_be(&data[sound_offset]);
               
               //adrs
               unsigned int adrs_offset = read_u32_be(&data[sound_offset+4]);
               if(adrs_offset != 0) {
                  sound_banks.banks[i].sounds[j].adrs = malloc(8 * sizeof(unsigned));
                  for (k = 0; k < 8; k++) {
                     sound_banks.banks[i].sounds[j].adrs[k] = read_u16_be(&data[adrs_offset+sound_bank_offset+16+k*2]);
                  }
               }
               else {
                  sound_banks.banks[i].sounds[j].adrs = NULL;
               }
               
               //wav_prev
               unsigned int wav_prev_offset = read_u32_be(&data[sound_offset+8]);
               if(wav_prev_offset != 0) {
                  sound_banks.banks[i].sounds[j].wav_prev = read_wave_table(data, wav_prev_offset + sound_bank_offset + 16, sound_bank_offset);
               }
               else {
                  sound_banks.banks[i].sounds[j].wav_prev = NULL;
               }
               unsigned int flt = read_u32_be(&data[sound_offset+12]);
               sound_banks.banks[i].sounds[j].key_base_prev = *((float*)&flt);
               
               //wav
               unsigned int wav_offset = read_u32_be(&data[sound_offset+16]);
               if(wav_offset != 0) {
                  sound_banks.banks[i].sounds[j].wav = read_wave_table(data, wav_offset + sound_bank_offset + 16, sound_bank_offset);
               }
               else {
                  sound_banks.banks[i].sounds[j].wav = NULL;
               }
               flt = read_u32_be(&data[sound_offset+20]);
               sound_banks.banks[i].sounds[j].key_base = *((float*)&flt);
               
               //wav_sec
               unsigned int wav_sec_offset = read_u32_be(&data[sound_offset+24]);
               if(wav_sec_offset != 0) {
                 sound_banks.banks[i].sounds[j].wav_sec = read_wave_table(data, wav_sec_offset + sound_bank_offset + 16, sound_bank_offset);
               }
               else {
                  sound_banks.banks[i].sounds[j].wav_sec = NULL;
               }
               flt = read_u32_be(&data[sound_offset+28]);
               sound_banks.banks[i].sounds[j].key_base_sec = *((float*)&flt);
            }
            else {
               sound_banks.banks[i].sounds[j].adrs = NULL;
               sound_banks.banks[i].sounds[j].wav_prev = NULL;
               sound_banks.banks[i].sounds[j].wav = NULL;
               sound_banks.banks[i].sounds[j].wav_sec = NULL;
            }
         }
       }
       
       //percussion
       if (sound_banks.banks[i].percussion_count > 0) {
         unsigned int perc_table_offset = read_u32_be(&data[sound_bank_offset+16]) + sound_bank_offset + 16;
         
          sound_banks.banks[i].percussions.items = malloc(sound_banks.banks[i].percussion_count * sizeof(percussion));
         for (j = 0; j < sound_banks.banks[i].percussion_count; j++) {
            unsigned int perc_offset = read_u32_be(&data[perc_table_offset+j*4]);
            
            if(perc_offset != 0)
            {
              perc_offset += sound_bank_offset + 16;
               sound_banks.banks[i].percussions.items[j].unknown_1 = data[perc_offset];
               sound_banks.banks[i].percussions.items[j].pan = data[perc_offset+1];
               sound_banks.banks[i].percussions.items[j].unknown_2 = read_u16_be(&data[perc_offset+2]);
               
               //wav
               unsigned int wav_offset = read_u32_be(&data[perc_offset+4]);
               if(wav_offset != 0) {
                 sound_banks.banks[i].percussions.items[j].wav = read_wave_table(data, wav_offset + sound_bank_offset + 16, sound_bank_offset);
               }
               unsigned int flt = read_u32_be(&data[perc_offset+8]);
               sound_banks.banks[i].percussions.items[j].key_base = *((float*)&flt);
               
               //adrs
               unsigned int adrs_offset = read_u32_be(&data[perc_offset+12]);
               if(adrs_offset != 0) {
                  sound_banks.banks[i].percussions.items[j].adrs = malloc(8 * sizeof(unsigned));
                 for (k = 0; k < 8; k++) {
                    sound_banks.banks[i].percussions.items[j].adrs[k] = read_u16_be(&data[adrs_offset+sound_bank_offset+16+k*2]);
                 }
               }
               else {
                  sound_banks.banks[i].percussions.items[j].adrs = NULL;
               }
            }
            else
            {
               sound_banks.banks[i].percussions.items[j].adrs = NULL;
               sound_banks.banks[i].percussions.items[j].wav = NULL;
            }
         }
       }
      }
   }
   
   return sound_banks;
}
   



// sfx standalone executable
#ifdef SFX_STANDALONE

int sfx_encode_file(const char *in_file, const char *out_file)
{
   FILE *in;
   FILE *out;
   unsigned char *in_buf = NULL;
   unsigned char *out_buf = NULL;
   size_t file_size;
   size_t bytes_read;
   int bytes_encoded;
   int bytes_written;
   int ret_val = 0;

   in = fopen(in_file, "rb");
   if (in == NULL) {
      return 1;
   }

   // allocate buffer to read entire contents of files
   fseek(in, 0, SEEK_END);
   file_size = ftell(in);
   fseek(in, 0, SEEK_SET);
   in_buf = malloc(file_size);

   // read bytes
   bytes_read = fread(in_buf, 1, file_size, in);
   if (bytes_read != file_size) {
      ret_val = 2;
   }
   else
   {
      
      //OKAY, WE NEED TO HAVE A .S FILE TO EXPORT THE CORRECT INFO INTO. WE ALSO NEED TO HAVE A 
      
      // allocate worst case length
      out_buf = malloc(MIO0_HEADER_LENGTH + ((file_size+7)/8) + file_size);

      // compress data in MIO0 format
      bytes_encoded = mio0_encode(in_buf, file_size, out_buf);

      // open output file
      out = fopen(out_file, "wb");
      if (out == NULL) {
         ret_val = 4;
         goto free_all;
      }
      else
      {

         // write data to file
         bytes_written = fwrite(out_buf, 1, bytes_encoded, out);
         if (bytes_written != bytes_encoded) {
            ret_val = 5;
         }

         // clean up
         fclose(out);
      }
   }

   if (out_buf) {
      free(out_buf);
   }
   if (in_buf) {
      free(in_buf);
   }
   fclose(in);

   return ret_val;
}




typedef struct
{
   char *in_filename;
   char *out_filename;
   unsigned long sampling_rate;
} arg_config;

static arg_config default_config =
{
   NULL,
   NULL
};

static void print_usage(void)
{
   ERROR("Usage: sfx-encode [-s SAMPLINGRATE] FILE [OUTPUT]\n"
         "\n"
         "sfx-encode: N64 SFX compression tool\n"
         "\n"
         "File arguments:\n"
         " -s SAMAPLINGRATE   sampling rate for sound data (default: 16000)\n"
         " FILE               input file\n"
         " [OUTPUT]           output file (default: FILE.out)\n");
   exit(1);
}

// parse command line arguments
static void parse_arguments(int argc, char *argv[], arg_config *config)
{
   int i;
   int file_count = 0;
   if (argc < 2) {
      print_usage();
      exit(1);
   }
   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
         switch (argv[i][1]) {
            case 's':
               config->sampling_rate = strtoul(argv[i], NULL, 0);
               break;
            default:
               print_usage();
               break;
         }
      } else {
         switch (file_count) {
            case 0:
               config->in_filename = argv[i];
               break;
            case 1:
               config->out_filename = argv[i];
               break;
            default: // too many
               print_usage();
               break;
         }
         file_count++;
      }
   }
   if (file_count < 1) {
      print_usage();
   }
}

int main(int argc, char *argv[])
{
   char out_filename[FILENAME_MAX];
   arg_config config;
   int ret_val;

   // get configuration from arguments
   config = default_config;
   parse_arguments(argc, argv, &config);
   if (config.out_filename == NULL) {
      config.out_filename = out_filename;
      sprintf(config.out_filename, "%s.out", config.in_filename);
   }

   // operation
   ret_val = sfx_encode_file(config.in_filename, config.out_filename);

   switch (ret_val) {
      /*case 1:
         ERROR("Error opening input file \"%s\"\n", config.in_filename);
         break;
      case 2:
         ERROR("Error reading from input file \"%s\"\n", config.in_filename);
         break;
      case 3:
         ERROR("Error decoding MIO0 data. Wrong offset (0x%X)?\n", config.offset);
         break;
      case 4:
         ERROR("Error opening output file \"%s\"\n", config.out_filename);
         break;
      case 5:
         ERROR("Error writing bytes to output file \"%s\"\n", config.out_filename);
         break;*/
   }

   return ret_val;
}
#endif // SFX_STANDALONE

