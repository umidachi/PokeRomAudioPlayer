#include <iostream>
#include <vector>
#include <cstdint>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <portaudio.h>

//#define SAMPLE_RATE 44100
//#define FRAMES_PER_BUFFER 64
#define SAMPLE_RATE 22100
#define FRAMES_PER_BUFFER 32


//const int SAMPLE_RATE = 44100;
double volume_rate =0.5;
double vibratoDepth = 0.0;
double vibratoRate = 0.0;
int noteSpeed = 0; // グローバルに持つ これを変える。
int tempo_value;
int bpm;
int mpqn;
double beatSec;
//int octave;
std::vector<int16_t> song;


struct Note { 
        double freq; 
        int beats; 
        int octave;
};
struct music_pointer{
    uint8_t bank;
    uint8_t channelNum;
    std::vector<size_t> addrList;

};

// ====== 状態構造体 ======
struct PlaybackState {
    const std::vector<uint8_t>* program; // ROMデータ
    size_t pc = 0;                       // プログラムカウンタ
    uint8_t bank;                        //プログラムカウンタ0の時のロムバンク
    size_t addr;                         //プログラムカウンタ0の時の相対アドレス
    uint8_t channel;                 // 仮: チャンネル番号
    uint8_t channelNum=1;               //チャネルのトータル数
    std::vector<int16_t> channelSong;    // 直近のノート波形
    size_t sampleIndex = 0;              // channelSong 再生中の位置
    uint8_t noteSpeed=1; // グローバルに持つ これを変える。

    uint8_t octave = 4;  
    uint8_t vibratorate = 0; // ビブラートの遅延
    uint8_t vibratoDepth = 0; // ビブラートの深さ
    float_t volume_rate=0;

    // サウンドコール用のスタック
    std::vector<size_t> callStack;
};


// ===== 音楽解析関数（簡略化） =====
double ticksToLen(int tks,int noteSpeed) {
    return (tks / 4.0) * (noteSpeed / 12.0) * beatSec;
}


std::vector<int16_t> generateTone(double volume_rate, double freq, double durationSec,
                                  double vibratoDepth = 0.0, double vibratoRate = 5.0) {
    std::vector<int16_t> samples;
    int totalSamples = static_cast<int>(SAMPLE_RATE * durationSec);
    samples.reserve(totalSamples);

    for (int i = 0; i < totalSamples; i++) {
        double t = (double)i / SAMPLE_RATE;
        double curFreq = freq;

        // Haskell版風のvibrato
        if (vibratoDepth > 0.0) {
            curFreq *= (1.0 + vibratoDepth * sin(2.0 * M_PI * vibratoRate * t));
        }

        double val = (fmod(t * curFreq, 1.0) < 0.5) ? volume_rate : -volume_rate;
        samples.push_back(static_cast<int16_t>(val * 32767));
    }

    return samples;
}
#include <vector>
#include <cstdint>
#include <cmath>
#include <random>

// ===================================================
// ゲームボーイ風ノイズ波生成
// ===================================================

std::vector<int16_t> generateDrumTone(double volume_rate, double durationSec, int drumType = 0) {
    std::vector<int16_t> samples;
    int totalSamples = static_cast<int>(SAMPLE_RATE * durationSec);
    samples.reserve(totalSamples);

    // ----------------------------
    // ドラムタイプ別パラメータ設定
    // ----------------------------
    // drumType: 0=Kick, 1=Snare, 2=HiHat
    int shiftClockFreq, dividingRatio, widthMode;
    double decayRate;

    switch (drumType) {
        case 0: // Kick
            shiftClockFreq = 10;
            dividingRatio = 3;
            widthMode = 0; // 15bit LFSR
            decayRate = 0.995;
            break;
        case 1: // Snare
            shiftClockFreq = 7;
            dividingRatio = 1;
            widthMode = 1; // 7bit LFSR
            decayRate = 0.992;
            break;
        case 2: // HiHat
        default:
            shiftClockFreq = 4;
            dividingRatio = 0;
            widthMode = 1;
            decayRate = 0.990;
            break;
    }

    // ----------------------------
    // LFSRノイズ生成
    // ----------------------------
    uint16_t lfsr = 0x7FFF; // 初期状態
    double amplitude = volume_rate * 32767.0;
    double env = 1.0;
    int counter = 0;

    // ゲームボーイ風のノイズ周波数計算
    // f = 524288 Hz / dividingRatio / 2^(shiftClockFreq+1)
    // dividingRatio==0 の場合は 0.5 として扱う
    double ratio = (dividingRatio == 0) ? 0.5 : dividingRatio;
    double noiseFreq = 524288.0 / ratio / std::pow(2.0, shiftClockFreq + 1);
    double samplesPerToggle = SAMPLE_RATE / noiseFreq;

    for (int i = 0; i < totalSamples; ++i) {
        if (counter++ >= samplesPerToggle) {
            counter = 0;
            // LFSR 更新
            uint16_t bit = (lfsr ^ (lfsr >> 1)) & 1;
            lfsr = (lfsr >> 1) | (bit << 14);
            if (widthMode)
                lfsr = (lfsr & ~0x40) | (bit << 6);
        }

        double value = (lfsr & 1) ? amplitude * env : -amplitude * env;
        samples.push_back(static_cast<int16_t>(value));

        env *= decayRate; // 短時間で音を減衰
        if (env < 0.001) break; // 無音なら終了
    }

    return samples;
}


void setTempo(const std::vector<uint8_t>& bytes, size_t pc) {
    //tempo ED xx xx = Set tempo
    uint16_t addr = bytes[pc+1] | (bytes[pc+2] << 8);

    tempo_value= bytes[pc+2]; //おそらくの計算式 1バイト目どこ行った・・・。
    //bpm = 14400 / tempo_value;
    
    if (tempo_value >= 100){
    mpqn = 60000000 / ( (362 - tempo_value) / 2 );
    }else if(tempo_value < 100){
    mpqn =60000000 / ((392 - tempo_value) / 2);
    }
    beatSec = mpqn / 1000000.0;

    std::cout << "--Tempo:" << std::dec << tempo_value<< std::endl;
    std::cout << "--set beatSec:"  << beatSec << std::endl;

}
void setVolume(const std::vector<uint8_t>& bytes, size_t pc) {
    
    uint8_t value = bytes[pc+1];
    std::cout << "set Volume to :" << std::hex << (int)value << std::endl;


}

void setDuty(const std::vector<uint8_t>& bytes, size_t pc) {
    uint8_t value = bytes[pc+1];
    std::cout << "set Duty to 0x" << std::hex << (int)value << std::endl;
    std::cout << "--DutyCycle:" << std::dec << (int)value << std::endl;

}
void setVibrato(const std::vector<uint8_t>& bytes, size_t pc ,PlaybackState& state) {
  //EA xx yz = Vibrato (X = Delay, Y = Depth, Z = Rate)
    uint8_t value1= bytes[pc+1]; //
    uint8_t Delay= value1; //

    uint8_t value2 = bytes[pc+2];
    uint8_t Byte = value2;

    state.vibratoDepth=((Byte >> 4) & 0x0F) * (1/127);    
    state.vibratorate = Byte & 0x0F * (1/127);

    std::cout << "set Vibrato to 0x" << std::hex << (int)value1 << std::hex << (int)value2 << std::endl;
    std::cout << "--Delay:" << std::dec << (int)Delay<< std::endl;
    std::cout << "--Depth:" << std::dec << (int)state.vibratoDepth<< std::endl;
    std::cout << "--Rate:" << std::dec << (int)state.vibratorate<< std::endl;    
}
void TogglePerfectPitch(const std::vector<uint8_t>& bytes, size_t pc) {
    uint8_t value1 = bytes[pc];
    std::cout << "Toggle Perfect Pitch to 0x" << std::hex << (int)value1 << std::endl;
}
void setNoteType(const std::vector<uint8_t>& bytes, size_t pc,PlaybackState& state) {
 
    uint8_t byte1 = bytes[pc];
    uint8_t byte2 = bytes[pc+1];
    
    // 上位4ビット = speed, 下位4ビット =
    state.noteSpeed = byte1 & 0x0F;
    int vol = (byte2 >> 4) & 0x0F;
    int fade = byte2 & 0x0F;

    state.volume_rate=vol/127.0;

    std::cout << "< CH" << std::hex << (int)state.channel  << " >" << std::endl;
    std::cout << "set NoteType 0x:" << std::hex << (int)byte1 << std::hex << (int)byte2 << std::endl;
    std::cout << "--Speed:" << std::dec << state.noteSpeed<< std::endl;
    std::cout << "--Vol:" << std::dec << vol<< std::endl;
    std::cout << "--Fade:" << std::dec << fade<< std::endl;


}
void setDrumSpeed(const std::vector<uint8_t>& bytes, size_t pc) {
  //EA xx yz = Vibrato (X = Delay, Y = Depth, Z = Rate)
  uint8_t value = bytes[pc];
  uint8_t value1 = bytes[pc+1];
    std::cout << "set DrumSpeed to 0x" << std::hex << (int)value << std::hex << (int)value1 << std::endl;
}
void setOctave(const std::vector<uint8_t>& bytes, size_t pc,PlaybackState& state) {
   uint8_t value = bytes[pc];
   uint8_t low  = value & 0x0F;        // 下位4ビット（オペランドなど）
 
   std::cout << "< CH" << std::hex << (int)state.channel  << " >" << std::endl;
   std::cout << "SetOctave to 0x" << std::hex << (int)value  << std::endl;
   std::cout << "--Channel " << (int)state.channel << " SetOctave to " << state.octave << std::endl;

    //octave = low;
	switch (low){
						case 0:
							state.octave = 8;
							break;
						case 1:
							state.octave = 7;
							break;
						case 2:
							state.octave = 6;
							break;
						case 3:
							state.octave = 5;
							break;
						case 4:
							state.octave = 4;
							break;
						case 5:
							state.octave = 3;
							break;
						case 6:
							state.octave = 2;
							break;
						case 7:
							state.octave = 1;
							break;
						}

}


void setRest(const std::vector<uint8_t>& bytes, size_t pc,std::vector<int16_t>&ch_song,PlaybackState &state) {
    uint8_t value = bytes[pc];
    uint8_t low  = value & 0x0F;        // 下位4ビット（オペランドなど）
    //double beatlen=low;
    double tks = low;
    double duration = ticksToLen(tks,state.noteSpeed);       

    int silenceSamples = (int)(SAMPLE_RATE * duration);
    ch_song.insert(ch_song.end(), silenceSamples, 0);

    std::cout << "< CH" << std::hex << (int)state.channel  << " >" << std::endl;
    std::cout << "SetRest to 0x" << std::hex << (int)value  << std::endl;
}



void setNote(const std::vector<uint8_t>& bytes, size_t pc,std::vector<int16_t>&ch_song,PlaybackState& state) {

    uint8_t value = bytes[pc];
    uint8_t high = (value >> 4) & 0x0F; // 上位4ビット
    uint8_t low  = (value & 0x0F) +1;        // 下位4ビット（オペランドなど）

    struct Note_type {
    const char* name;   // ドレミ名
        double frequency;   // 周波数(Hz)
    };
    

    uint8_t noteIndex=high;
    double beatlen=low;

    std::cout << "< CH" << std::hex << (int)state.channel  << " >" << std::endl;
    std::cout << "SetNote to 0x" << std::hex << (int)value  << std::endl;
     double tks = beatlen ;
    double duration = ticksToLen(tks, state.noteSpeed);       
            int ch=1;   //暫定的な設定
    
             // 基準オクターブ4（noteTableの周波数はC4〜B4想定）との相対オクターブ差を計算
             // haskel版を採用
            // チャンネル補正: Ch3は+1、それ以外は+2
            //int channelOffset = (ch == 3) ? 1 : 2;
             int midiNote = noteIndex + (12 * (state.octave + 2 ))  ;
            double freqWithOct = 440.0 * pow(2.0, (midiNote - 69) / 12.0);
 
            auto tone = generateTone(state.volume_rate, freqWithOct, duration, state.vibratoDepth, state.vibratorate); // vibratoあり
           //auto tone = generateTone(state.volume_rate, freqWithOct, duration, 0, 0); // vibratoなし
            std::cout << "--calFreq = " << freqWithOct << " Hz" << std::endl;      
            std::cout << "--duration = " << duration << "s" << std::endl;      


            ch_song.insert(ch_song.end(), tone.begin(), tone.end());
       

}

int setSoundloop(const std::vector<uint8_t>& bytes, size_t pc,PlaybackState& state) {
    //ex: fe 00 30 7a
    // 00は無限ループ 30 7aがリトルエンディアンのアドレス
    uint8_t low = bytes[pc+2];
    uint8_t high = bytes[pc+3];
    uint16_t jump_label = low |(high << 8);

    uint16_t calcpc= jump_label-state.addr;   //ジャンプラベルのアドレスと、スタートの位置のアドレス の差分を算出し、ジャンプ先のPCを計算

    std::cout << "< CH" << std::hex << (int)state.channel  << " >" << std::endl;
    std::cout << "SetSoundloop to 0x" << std::hex << (int)jump_label  << std::endl;


    return calcpc;
}

// 指定した開始位置から特定のバイトパターンが出るまで抽出
std::vector<uint8_t> extractSection(
    const std::vector<uint8_t>& rom,       // ROM全体
    size_t startOffset,                    // 開始アドレス
    const std::vector<uint8_t>& pattern    // 終了パターン
) {
    if (startOffset >= rom.size()) {
        throw std::runtime_error("startOffset が ROMサイズを超えています");
    }

    std::vector<uint8_t> section;
    for (size_t i = startOffset; i < rom.size(); ++i) {
        section.push_back(rom[i]);

        // 終了パターンチェック
        if (i + pattern.size() <= rom.size()) {
            bool match = true;
            for (size_t j = 0; j < pattern.size(); ++j) {
                if (rom[i + j] != pattern[j]) {
                    match = false;
                    break;
                }
            }
            if (match) {
                // パターンを含めて終了
                section.insert(section.end(), pattern.begin(), pattern.end());
                break;
            }
        }
    }

    return section;
}

size_t setSoundCall(const std::vector<uint8_t>& bytes, size_t pc, PlaybackState& state) {
    // FD xx xx = サウンドコール
    uint8_t low  = bytes[pc+1];
    uint8_t high = bytes[pc+2];
    uint16_t callAddr = low | (high << 8);

    // 呼び出し元をスタックに保存（FDの次の命令の位置）
    state.callStack.push_back(pc + 3);

    // 実際のROM内オフセットを計算
    size_t targetPc = callAddr - state.addr;

    std::cout << "< CH" << std::hex << (int)state.channel  << " >" << std::endl;
    std::cout << "SoundCall to 0x" << std::hex << callAddr 
              << " (jump PC=" << targetPc << ")" << std::endl;

    // ジャンプ先のPCを返す
    state.channelSong.clear();
    return targetPc;
}

size_t decodeByteNote( PlaybackState& state) 
{
    size_t pc = state.pc;
    state.channelSong.clear();   // 呼び出すたびに初期化
    std::vector<uint8_t> program = *state.program;

    while (1) {
        uint8_t byte = program[pc];

        if (byte == 0xED) { // tempo
            setTempo(program, pc);
            pc += 3;
        }
        else if (byte == 0xF0) { // volume
            setVolume(program, pc);
            pc += 2;
        }
        else if (byte == 0xEC) { // duty
            setDuty(program, pc);
            pc += 2;
        }
        else if (byte == 0xEA) { // vibrato
            setVibrato(program, pc,state);
            pc += 3;
        }
        else if (byte == 0xE8) { // perfect pitch
            TogglePerfectPitch(program, pc);
            pc += 1;
        }else if(byte ==0xFE){
            pc=setSoundloop(program,pc,state);
            //return pc;
        }else if(byte == 0xFD){
            pc=setSoundCall(program,pc,state);
            //return pc;
        }else if(byte == 0xFF){ //End of song or return from macro          
           if (!state.callStack.empty()) {
                pc = state.callStack.back();  // 戻り先を取り出す
                state.callStack.pop_back();         // 使い終わったのでスタックから削除
            }else{
                pc=0;
            }
           

        }
        
        else {
            uint8_t high = (byte >> 4) & 0x0F;

            switch (high) {
                case 0xC: // Rest
                    setRest(program, pc, state.channelSong, state);
                    pc += 1;
                    return pc;

                case 0xD: // Note type / Drum speed
                    if (state.channel != 4) {
                        setNoteType(program, pc,state);
                        pc += 1;
                    } else {
                        setDrumSpeed(program, pc);
                    }
                    break;

                case 0xE: // Octave
                    setOctave(program, pc,state);
                    break;

                case 0xF: // Special command
                    //pc += 1;
                    break;

                default:  // Note
                    setNote(program, pc, state.channelSong,state);
                    pc += 1;
                    return pc;
            }
            pc++;
        }
    }
}

// ===== PortAudio コールバック =====
static int paCallback(const void* inputBuffer, void* outputBuffer,
                      unsigned long framesPerBuffer,
                      const PaStreamCallbackTimeInfo* timeInfo,
                      PaStreamCallbackFlags statusFlags,
                      void* userData) {
    auto* states = reinterpret_cast<PlaybackState*>(userData);
    PlaybackState& ch1 = states[0];
    PlaybackState& ch2 = states[1];
    PlaybackState& ch3 = states[2];
    PlaybackState& ch4 = states[3];

    std::vector<PlaybackState*> channels = {&ch1, &ch2, &ch3, &ch4};
    //std::vector<PlaybackState*> channels = {&ch1, &ch2, &ch3};
    //uint8_t max_ch=channels[0]->channelNum;
    uint8_t max_ch=3;   
    float* out = reinterpret_cast<float*>(outputBuffer);
    float sample[] = {0,0,0,0};
    //float sample[] = {0,0,0};
    float mixed=0;
   
        for (unsigned int i = 0; i < framesPerBuffer; i++) {
            for(int n=0;n<max_ch;n++){
            PlaybackState* ch = channels[n];
            if (ch->channelSong.empty() && ch->pc < ch->program->size()){
                ch->pc = decodeByteNote(*ch);
            }
            sample[n]= 0.0f;
            if (!ch->channelSong.empty()) {
                if (ch->sampleIndex < ch->channelSong.size()) {
                    sample[n]= static_cast<float>(ch->channelSong[ch->sampleIndex]) / 32768.0f;
                    ch->sampleIndex++;
                } else {
                    ch->channelSong.clear();
                    ch->sampleIndex = 0;
                }
            }
        //mixed = (sample[0]+sample[1]+sample[2]+sample[3]) /3.3f;
            mixed = (sample[0]+sample[1]+sample[2]+sample[3]) /max_ch;


        }
        // 合成（単純平均でもよい）

        //float mixed = sample4;
        
        // 左右に出力
        //*out++ = sample1;
        //*out++ = sample2;
        *out++ = mixed;
        *out++ = mixed;
    }
    // 両方終わったら終了
/*    if (ch1.pc >= ch1.program->size() && ch1.channelSong.empty() &&
        ch2.pc >= ch2.program->size() && ch2.channelSong.empty()&&
        ch3.pc >= ch3.program->size() && ch3.channelSong.empty()&&
        ch4.pc >= ch4.program->size() && ch4.channelSong.empty()) {
        return paComplete;
    }
*/      

    return paContinue;
}


// ===== バイナリ読み込み =====
std::vector<uint8_t> readBinaryFile(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) throw std::runtime_error("ファイルを開けません");

    file.seekg(0, std::ios::end);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<uint8_t> buffer(size);
    if (!file.read(reinterpret_cast<char*>(buffer.data()), size))
        throw std::runtime_error("読み込み失敗");

    return buffer;
}

void getMusicPointer(const std::vector<uint8_t>* rom,size_t bank, size_t pointer_address, std::vector<music_pointer>& musicPointerList){
    //ポインタのアドレスから、音楽データのアドレスを取得する。
    std::cout << "bank:0x" << std::hex << (int)bank << std::endl;
    size_t targetAddress = (bank << 14) + (pointer_address & 0x3fff);
    std::cout << "targetAddress:0x" << std::hex << targetAddress<< std::endl;

    for(int i=0;i<3;i++){

       int Header= (int)(*rom)[targetAddress+i];
        std::cout << (int)Header<< std::endl;

        if(Header != 0xff){
            std::cout << "no sfx header: FF FF FF"<< std::endl;
            return;
        }
    }
    std::cout <<"get Sfx header" << std::endl;
    //音楽データのアドレスであるとあることを確認したので、アドレスを取得する。
    int i=0;
    int pc = i+3;//FF FF FF の次の1バイトを読み取る
    int endpoint;

    if(bank==0x1f){
        endpoint=0x7FFF-0x4000; //bank1fの終端アドレス
    }else if(bank==0x02){
        endpoint=0xBFFF-0x4000; //bank02の終端アドレス
    }else{
        endpoint=500;
    }

    while(pc < endpoint){

        uint8_t byte  = (*rom)[targetAddress + pc];

        uint8_t high = (byte >> 4) & 0x0F;
        uint8_t low  = (byte & 0x0F);       

        int channelNum;


        switch(high){
            case 0:
                channelNum=1;
                break;
            case 4:
                channelNum=2;    
                break;
            case 8:
                channelNum=3;
                break;
            case 0xC:
                channelNum=4;
                break;
            default:
                std::cout << "Unknown Channel Num" << std::endl;
                return;
        }
        music_pointer mp;
        mp.bank = bank;
        mp.channelNum=channelNum;
        pc=pc+1;
        
        //次のバイトから最初に呼んだチャネル数分のポインタを取得する。

        for(int i=0;i<channelNum;i++){
            uint8_t low_byte;
            uint8_t high_byte; 
            uint8_t header_byte;
            //2バイトを読み取り、アドレスに変換する
            if(i==0){
                low_byte = (*rom)[targetAddress + pc];
                high_byte  = (*rom)[targetAddress + (pc + 1)];    
                pc=pc+2;
            }else{
                header_byte = (*rom)[targetAddress + pc];
                low_byte = (*rom)[targetAddress + pc + 1];
                high_byte  = (*rom)[targetAddress + pc + 2];
                pc=pc+3;
            }
            
            uint16_t addr = low_byte | (high_byte << 8);
            //size_t startOffset1 = (0x1F << 14) + (addr & 0x3fff);
            mp.addrList.push_back(addr);

        }
        /*    std::cout << "------------------------" << std::endl;
            std::cout << "bank: " << std::hex << (int)mp.bank << std::endl;
            std::cout << "ch_num: " << std::hex << (int)mp.channelNum << std::endl;
            for(int i=0;i<mp.channelNum;i++){
            std::cout << "ch"<< i << ":"  <<std::hex << mp.addrList[i] << std::endl;
            }
        */
            musicPointerList.push_back(mp);
    }

    std::cout << "Total Music Num:" << musicPointerList.size() << std::endl;

}


// ===== main =====
int main(int argc, char* argv[]) {
 if (argc < 2) {
        std::cerr << "使い方: " << argv[0] << " 音楽リスト番号 " << std::endl;
        return 1;
    }

    int num = std::atof(argv[1]);

    try {
        std::vector<uint8_t> rom = readBinaryFile("./pokeblue.gbc");
        std::vector<music_pointer> musicPointerList;


        getMusicPointer(&rom,0x1f,0x4000,musicPointerList);
        getMusicPointer(&rom,0x02,0x4000,musicPointerList);
        getMusicPointer(&rom,0x08,0x4000,musicPointerList);

        std::cout << "Total Music Num:" << musicPointerList.size() << std::endl;

        // ★ 再生状態を用意
        PlaybackState state_ch1, state_ch2, state_ch3, state_ch4;

        std::vector<PlaybackState>state_ch={state_ch1,state_ch2,state_ch3,state_ch4};
        std::vector<uint8_t> list = {0xe6,0xe3,0x6c, 0xeb, 0xd9,0xd8,0xda,0x6f,0x6d,0x70,0x71,0x72,0x73,0x74,0x75,0x76,0x77}; // 終了パターン
        uint8_t  N=list[num]; // 再生する音楽の番号
        std::vector<uint8_t> program[4];

        for(int i=0;i<musicPointerList[N].channelNum;i++){
            int Offset;

            program[i].clear();
            state_ch[i].bank=musicPointerList[N].bank;
            state_ch[i].addr=musicPointerList[N].addrList[i];
            state_ch[i].channelNum=musicPointerList[N].channelNum;
            Offset= (state_ch[i].bank << 14) + (state_ch[i].addr & 0x3fff);  
            program[i] = {rom.begin() + Offset, rom.begin() + Offset + 500}; //300とかだと止まる。
            state_ch[i].program = &program[i];
            state_ch[i].channel=i+1;
    }

        //デバッグ用　曲番号がわからないときに使う。
        for(int i=0;i<musicPointerList.size();i++){
            std::cout << "Music bank:" << (int)musicPointerList[i].bank << std::endl;
            std::cout << "Music Num:" << i << std::endl;
            std::cout << "ch_num: " << std::hex << (int)musicPointerList[i].channelNum << std::endl;
            for(int j=0;j<musicPointerList[i].channelNum;j++){
                std::cout << "ch"<< j << ":"  <<std::hex << musicPointerList[i].addrList[j] << std::endl;
            }
        }



        Pa_Initialize();
        PaStream* stream;
        PlaybackState states[4] = {state_ch[0], state_ch[1],state_ch[2],state_ch[3]};  


        Pa_OpenDefaultStream(&stream, 0, 2, paFloat32, SAMPLE_RATE,
                             FRAMES_PER_BUFFER, paCallback, &states);
        Pa_StartStream(stream);

        std::cout << "Playing... Enter to Stop" << std::endl;
        std::cin.get();

        Pa_StopStream(stream);
        Pa_CloseStream(stream);
        Pa_Terminate();
    }
    catch (const std::exception& e) {
        std::cerr << "エラー: " << e.what() << std::endl;
        return 1;
    }
}
